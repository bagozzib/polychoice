// firth.cpp -- Firth's penalty for mnl() (mode 0) and rologit() (mode 1): the
// information I(theta) of the Jeffreys-prior penalty 1/2 log|I| and, given
// I^{-1}, each case's share of the penalty's gradient.
//
// A stage chooses from a remaining set R with the utilities of parameter block
// b (rank_test()'s split model has two blocks). The multinomial / conditional
// logit has one stage per case, R its choice set; a ranking has one per ranked
// position with an alternative below it (Beggs, Cardell and Hausman 1981).
//   Untied stage: D_R' W_R D_R, W_R = diag(P_R) - P_R P_R' (the information of
//   the stage given R; the observed information of the ranking model).
//   Tie group of d alternatives (exact likelihood, Allison and Christakis 1994,
//   eq. 8): the probability that the group is the top-d set of R; the group is
//   one categorical outcome over the size-d subsets S of R, and contributes its
//   expected information D' Omega D, Omega = sum_S pi_S g_S g_S',
//   g_S = d log pi_S / dV. pi_S and its first and second V-derivatives come from
//   the recursion over the subsets U of R with |U| <= d,
//     f(U) = sum_{a in U} f(U \ a) q_t(a | R \ (U \ a)),   t = |U|,
//   q_t the conditional-logit probability at the group's t-th stage (positive
//   terms only; at d = 1 Omega = W_R).
// Gradient of 1/2 log|I| (H = D I^{-1} D', h = diag(H)):
//   untied stage  1/2 D' W a,  a = h - 2 H P
//   tie group     1/2 D' sum_S pi_S [(g_S' H g_S) g_S + 2 G_S H g_S],
//                 G_S = d^2 log pi_S / dV dV'.
// Derivations and their numerical checks: DESIGN.md (section 1, Phase 4) and
// verification/probes/2026-11-05. Design rows are handled as sparse vectors,
// as in mnl_ll.cpp and rologit_ll.cpp.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <cstdint>
using namespace Rcpp;

namespace {

typedef std::vector< std::pair<int, double> > SRow;

struct FSpec {
  int n, J, p, q, base, K0, nb, split, K, mode;
  std::vector<int> bpos;
  int gpos;
};

FSpec make_fspec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W, const NumericMatrix& off,
                 const IntegerMatrix& avail, const IntegerMatrix& R, int mode, int J, int base, int split) {
  FSpec s;
  s.n = X.nrow(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.base = base; s.split = split; s.mode = mode;
  s.bpos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  s.K0 = pos; s.nb = (split > 0) ? 2 : 1; s.K = s.nb * s.K0;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (W.nrow() != s.n * J) stop("internal: W must have n * J rows.");
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  if (R.nrow() != s.n || R.ncol() != J) stop("internal: the rank matrix does not match the design.");
  return s;
}

inline int block_at(const FSpec& s, int position) { return (s.split > 0 && position > s.split) ? 1 : 0; }

void utilities(int i, const FSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
               const NumericMatrix& off, const IntegerMatrix& avail, std::vector<double>& V) {
  const int J = s.J;
  for (int b = 0; b < s.nb; ++b) {
    const int o = b * s.K0;
    for (int j = 0; j < J; ++j) {
      if (!avail(i, j)) { V[b * J + j] = R_NegInf; continue; }
      double v = off(i, j);
      if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) v += X(i, c) * theta[o + s.bpos[j] + c];
      for (int c = 0; c < s.q; ++c) v += W(i * J + j, c) * theta[o + s.gpos + c];
      V[b * J + j] = v;
    }
  }
}

void design_row(int i, int b, int j, const FSpec& s, const NumericMatrix& X, const NumericMatrix& W, SRow& out) {
  out.clear();
  const int o = b * s.K0;
  if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) { const double v = X(i, c); if (v != 0.0) out.push_back(std::make_pair(o + s.bpos[j] + c, v)); }
  for (int c = 0; c < s.q; ++c) { const double v = W(i * s.J + j, c); if (v != 0.0) out.push_back(std::make_pair(o + s.gpos + c, v)); }
}

// a stage group: the remaining set R (the group's alternatives and those below
// it), the group's size d, and the positions placed above it
struct Group { std::vector<int> R; int d; int pos; };

void case_groups(int i, const FSpec& s, const IntegerMatrix& avail, const IntegerMatrix& Rk, std::vector<Group>& out) {
  out.clear();
  const int J = s.J;
  if (s.mode == 0) {
    Group g; g.d = 1; g.pos = 0;
    for (int j = 0; j < J; ++j) if (avail(i, j)) g.R.push_back(j);
    if (g.R.size() >= 2) out.push_back(g);
    return;
  }
  std::vector< std::pair<int, int> > ranked;
  std::vector<int> unranked;
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    if (Rk(i, j) > 0) ranked.push_back(std::make_pair(Rk(i, j), j)); else unranked.push_back(j);
  }
  std::sort(ranked.begin(), ranked.end());
  int pos = 0;
  size_t a = 0;
  while (a < ranked.size()) {
    size_t e = a;
    while (e < ranked.size() && ranked[e].first == ranked[a].first) ++e;
    Group g; g.d = (int) (e - a); g.pos = pos;
    for (size_t t = a; t < e; ++t) g.R.push_back(ranked[t].second);
    const size_t nT = g.R.size();
    for (size_t t = e; t < ranked.size(); ++t) g.R.push_back(ranked[t].second);
    for (size_t t = 0; t < unranked.size(); ++t) g.R.push_back(unranked[t]);
    if (g.R.size() == nT) break;                     // nothing below: the group's term is 1, no stage
    out.push_back(g);
    pos += g.d;
    a = e;
  }
}

// the local coordinates of a group: (block index bi, alternative index a in R)
// -> bi * r + a, over the distinct blocks its stages use
struct Local {
  int r, nbu, m;
  std::vector<int> blocks;      // distinct blocks, in order of first use
  std::vector<int> stage_bi;    // block index of stage t = 1..d (entry t - 1)
};

Local make_local(const FSpec& s, const Group& g) {
  Local L; L.r = (int) g.R.size();
  for (int t = 1; t <= g.d; ++t) {
    const int b = block_at(s, g.pos + t);
    int bi = -1;
    for (size_t u = 0; u < L.blocks.size(); ++u) if (L.blocks[u] == b) bi = (int) u;
    if (bi < 0) { L.blocks.push_back(b); bi = (int) L.blocks.size() - 1; }
    L.stage_bi.push_back(bi);
  }
  L.nbu = (int) L.blocks.size(); L.m = L.nbu * L.r;
  return L;
}

// the conditional-logit probabilities of stage block b over the alternatives
// of R not in `taken` (P[a] = 0 for taken ones)
void stage_probs(const FSpec& s, const Group& g, const std::vector<double>& V, int b, uint64_t taken, std::vector<double>& P) {
  const int r = (int) g.R.size();
  // a mask covers at most 62 alternatives (tie_group refuses more); an untied
  // stage passes 0 and may have more than 64, where shifting would be
  // undefined behaviour (local review 2026-11-06)
  auto is_taken = [taken](int a) { return a < 64 && ((taken >> a) & 1ULL); };
  double vmax = R_NegInf;
  for (int a = 0; a < r; ++a) if (!is_taken(a)) vmax = std::max(vmax, V[b * s.J + g.R[a]]);
  double den = 0.0;
  for (int a = 0; a < r; ++a) { P[a] = is_taken(a) ? 0.0 : std::exp(V[b * s.J + g.R[a]] - vmax); den += P[a]; }
  for (int a = 0; a < r; ++a) P[a] /= den;
}

double n_subsets(int r, int d) {
  double tot = 0.0, c = 1.0;
  for (int t = 0; t <= d; ++t) { tot += c; c = c * (r - t) / (t + 1); }
  return tot;
}

// the largest number of size-t subsets of r elements over t = 1..d
double max_level(int r, int d) {
  double mx = 0.0, c = 1.0;
  for (int t = 0; t < d; ++t) { c = c * (r - t) / (t + 1); mx = std::max(mx, c); }
  return mx;
}

// A tie group: Omega (m x m, local coordinates) and, when H (m x m) is given,
// cvec = sum_S pi_S [(g'Hg) g + 2 G H g] (length m). The second derivatives F2
// (m x m per subset) enter cvec only, so they are computed only when H is
// given; their memory is bounded (2.5e7 doubles, 200 MB, per level) as the
// number of subsets is (1e5). `case_no` is the case's position in the
// estimation sample.
void tie_group(int case_no, const FSpec& s, const Group& g, const Local& L, const std::vector<double>& V,
               const double* H, std::vector<double>& Omega, std::vector<double>& cvec) {
  const int r = L.r, m = L.m, d = g.d;
  const bool need2 = (H != NULL);
  if (r > 62) stop("Under firth = TRUE or a prior, the %d-th case of the estimation sample has a tie group among more than 62 remaining alternatives; such a tie is too large for the penalty.", case_no);
  if (n_subsets(r, d) > 1e5)
    stop("Under firth = TRUE or a prior, the %d-th case of the estimation sample has a tie group of %d among %d remaining alternatives (more than 1e5 subsets to enumerate); this tie is too large for the penalty.", case_no, d, r);
  if (need2 && max_level(r, d) * (double) m * (double) m > 2.5e7)
    stop("Under firth = TRUE, the %d-th case of the estimation sample has a tie group of %d among %d remaining alternatives whose penalty gradient needs more than 200 MB; this tie is too large for the penalty.", case_no, d, r);
  Omega.assign((size_t) m * m, 0.0);
  cvec.assign(m, 0.0);
  // level t - 1: masks, index, f, F1 (m per mask), F2 (m * m per mask, with H only)
  const size_t m2 = need2 ? (size_t) m * m : 0;
  std::vector<uint64_t> pm(1, 0ULL);
  std::unordered_map<uint64_t, int> pidx; pidx[0ULL] = 0;
  std::vector<double> pf(1, 1.0), pF1((size_t) m, 0.0), pF2(m2, 0.0);
  std::vector<double> P(r), gl(m), F1p, F2p;
  for (int t = 1; t <= d; ++t) {
    const int b = L.blocks[L.stage_bi[t - 1]], o = L.stage_bi[t - 1] * r;
    std::vector<uint64_t> cm;
    for (uint64_t U = (1ULL << t) - 1ULL; U < (1ULL << r); ) {          // the size-t subsets (Gosper's hack)
      cm.push_back(U);
      const uint64_t c = U & (~U + 1ULL), rr = U + c;
      U = (((rr ^ U) >> 2) / c) | rr;
    }
    const size_t nc = cm.size();
    std::vector<double> cf(nc, 0.0), cF1(nc * m, 0.0), cF2(nc * m2, 0.0);
    std::unordered_map<uint64_t, int> cidx;
    for (size_t u = 0; u < nc; ++u) {
      if (u % 1000 == 999) Rcpp::checkUserInterrupt();
      const uint64_t U = cm[u];
      cidx[U] = (int) u;
      double* F1 = &cF1[u * m]; double* F2 = need2 ? &cF2[u * m2] : NULL;
      for (int a = 0; a < r; ++a) {
        if (!((U >> a) & 1ULL)) continue;
        const uint64_t prev = U & ~(1ULL << a);
        const int pi = pidx[prev];
        const double fp = pf[pi];
        if (fp <= 0.0) continue;
        stage_probs(s, g, V, b, prev, P);
        const double qa = P[a];
        if (qa <= 0.0) continue;
        std::fill(gl.begin(), gl.end(), 0.0);
        for (int c = 0; c < r; ++c) if (!((prev >> c) & 1ULL)) gl[o + c] = (c == a ? 1.0 : 0.0) - P[c];
        const double* F1q = &pF1[(size_t) pi * m];
        cf[u] += fp * qa;
        for (int x = 0; x < m; ++x) F1[x] += qa * (F1q[x] + fp * gl[x]);
        if (!need2) continue;
        const double* F2q = &pF2[(size_t) pi * m2];
        for (int x = 0; x < m; ++x) for (int y = 0; y < m; ++y)
          F2[x * m + y] += qa * (F2q[x * m + y] + F1q[x] * gl[y] + gl[x] * F1q[y] + fp * gl[x] * gl[y]);
        // the conditional logit's own curvature: -(diag P - P P') in block b
        for (int x = 0; x < r; ++x) {
          if ((prev >> x) & 1ULL) continue;
          for (int y = 0; y < r; ++y) {
            if ((prev >> y) & 1ULL) continue;
            F2[(o + x) * m + o + y] -= qa * fp * P[x] * ((x == y ? 1.0 : 0.0) - P[y]);
          }
        }
      }
    }
    pm.swap(cm); pidx.swap(cidx); pf.swap(cf); pF1.swap(cF1); pF2.swap(cF2);
  }
  // the size-d subsets: Omega and the gradient vector
  std::vector<double> gS(m), hg(m), Ghg(m);
  for (size_t u = 0; u < pm.size(); ++u) {
    const double f = pf[u];
    if (!(f > 1e-300)) continue;
    const double* F1 = &pF1[u * m]; const double* F2 = need2 ? &pF2[u * m2] : NULL;
    for (int x = 0; x < m; ++x) gS[x] = F1[x] / f;
    for (int x = 0; x < m; ++x) for (int y = 0; y < m; ++y) Omega[x * m + y] += f * gS[x] * gS[y];
    if (H) {
      double ghg = 0.0;
      for (int x = 0; x < m; ++x) { double v = 0.0; for (int y = 0; y < m; ++y) v += H[x * m + y] * gS[y]; hg[x] = v; ghg += gS[x] * v; }
      for (int x = 0; x < m; ++x) {                              // G hg, G = F2 / f - g g'
        double v = 0.0;
        for (int y = 0; y < m; ++y) v += (F2[x * m + y] / f - gS[x] * gS[y]) * hg[y];
        Ghg[x] = v;
      }
      for (int x = 0; x < m; ++x) cvec[x] += f * (ghg * gS[x] + 2.0 * Ghg[x]);
    }
  }
}

// the design rows of a group's local coordinates
void local_rows(int i, const FSpec& s, const Group& g, const Local& L, const NumericMatrix& X, const NumericMatrix& W, std::vector<SRow>& D) {
  D.assign(L.m, SRow());
  for (int bi = 0; bi < L.nbu; ++bi)
    for (int a = 0; a < L.r; ++a) design_row(i, L.blocks[bi], g.R[a], s, X, W, D[bi * L.r + a]);
}

} // namespace

// The information I(theta) of the penalty (K x K), weighted.
// [[Rcpp::export]]
NumericMatrix pc_firth_info_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                                IntegerMatrix R, int mode, int J, int base, int split, NumericVector wt) {
  FSpec s = make_fspec(theta, X, W, off, avail, R, mode, J, base, split);
  if (wt.size() != s.n) stop("internal: weights do not match the cases.");
  NumericMatrix I(s.K, s.K);
  std::vector<double> V(s.nb * J), P(J), Omega, cvec;
  std::vector<Group> groups;
  std::vector<SRow> D;
  for (int i = 0; i < s.n; ++i) {
    if (i % 1000 == 999) Rcpp::checkUserInterrupt();
    const double wi = wt[i];
    if (wi == 0.0) continue;
    utilities(i, s, theta, X, W, off, avail, V);
    case_groups(i, s, avail, R, groups);
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      const Group& g = groups[gi];
      Local L = make_local(s, g);
      const int m = L.m;
      if (g.d == 1) {
        P.assign(L.r, 0.0);
        stage_probs(s, g, V, L.blocks[0], 0ULL, P);
        Omega.assign((size_t) m * m, 0.0);
        for (int x = 0; x < m; ++x) for (int y = 0; y < m; ++y) Omega[x * m + y] = P[x] * ((x == y ? 1.0 : 0.0) - P[y]);
      } else tie_group(i + 1, s, g, L, V, NULL, Omega, cvec);
      local_rows(i, s, g, L, X, W, D);
      for (int x = 0; x < m; ++x) {
        if (D[x].empty()) continue;
        for (int y = 0; y < m; ++y) {
          if (D[y].empty()) continue;
          const double h = wi * Omega[x * m + y];
          if (h == 0.0) continue;
          for (size_t a = 0; a < D[x].size(); ++a)
            for (size_t b = 0; b < D[y].size(); ++b)
              I(D[x][a].first, D[y][b].first) += h * D[x][a].second * D[y][b].second;
        }
      }
    }
  }
  return I;
}

// Each case's share of the gradient of 1/2 log|I| (n x K, unweighted: the
// gradient is sum_i w_i row_i), given Iinv = I^{-1}.
// [[Rcpp::export]]
NumericMatrix pc_firth_adj_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                               IntegerMatrix R, int mode, int J, int base, int split, NumericMatrix Iinv) {
  FSpec s = make_fspec(theta, X, W, off, avail, R, mode, J, base, split);
  if (Iinv.nrow() != s.K || Iinv.ncol() != s.K) stop("internal: Iinv must be K x K.");
  NumericMatrix A(s.n, s.K);
  std::vector<double> V(s.nb * J), P, Omega, cvec, H, tu;
  std::vector<Group> groups;
  std::vector<SRow> D;
  for (int i = 0; i < s.n; ++i) {
    if (i % 1000 == 999) Rcpp::checkUserInterrupt();
    utilities(i, s, theta, X, W, off, avail, V);
    case_groups(i, s, avail, R, groups);
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      const Group& g = groups[gi];
      Local L = make_local(s, g);
      const int m = L.m;
      local_rows(i, s, g, L, X, W, D);
      // H = D I^{-1} D' over the local coordinates
      H.assign((size_t) m * m, 0.0);
      tu.assign(s.K, 0.0);
      for (int x = 0; x < m; ++x) {
        if (D[x].empty()) continue;
        std::fill(tu.begin(), tu.end(), 0.0);
        for (size_t a = 0; a < D[x].size(); ++a) {
          const int col = D[x][a].first; const double v = D[x][a].second;
          for (int k = 0; k < s.K; ++k) tu[k] += Iinv(k, col) * v;
        }
        for (int y = 0; y < m; ++y) {
          double hv = 0.0;
          for (size_t b = 0; b < D[y].size(); ++b) hv += tu[D[y][b].first] * D[y][b].second;
          H[x * m + y] = hv;
        }
      }
      if (g.d == 1) {
        P.assign(L.r, 0.0);
        stage_probs(s, g, V, L.blocks[0], 0ULL, P);
        // c = W a, a = h - 2 H P
        std::vector<double> av(m);
        for (int x = 0; x < m; ++x) { double hp = 0.0; for (int y = 0; y < m; ++y) hp += H[x * m + y] * P[y]; av[x] = H[x * m + x] - 2.0 * hp; }
        double pa = 0.0; for (int x = 0; x < m; ++x) pa += P[x] * av[x];
        cvec.assign(m, 0.0);
        for (int x = 0; x < m; ++x) cvec[x] = P[x] * (av[x] - pa);
      } else tie_group(i + 1, s, g, L, V, H.data(), Omega, cvec);
      for (int x = 0; x < m; ++x) {
        if (cvec[x] == 0.0) continue;
        for (size_t a = 0; a < D[x].size(); ++a) A(i, D[x][a].first) += 0.5 * cvec[x] * D[x][a].second;
      }
    }
  }
  return A;
}
