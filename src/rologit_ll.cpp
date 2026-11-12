// rologit_ll.cpp -- the likelihood engine of rologit(): the rank-ordered
// (exploded) logit for complete, partial, and tied rankings.
//
// Model (case i = 0..n-1, alternatives j = 0..J-1, base b, availability a_ij):
//   V_ij = off_ij + x_i' beta_j + w_ij' gamma,   beta_b = 0     (as in mnl_ll.cpp)
//   R_ij > 0: the rank of alternative j (smaller is preferred; equal values are
//   ties); R_ij = 0: available but unranked (tied below every ranked one).
//   The ranked alternatives form groups G_1, ..., G_K of equal rank. A group T
//   of d alternatives with B the alternatives ranked below it (the later groups
//   and the unranked) contributes (Allison and Christakis 1994, eq. 8)
//     sum over the orders p of T of
//       prod_{s=1}^{d} exp(V_{p_s}) / (sum_{t >= s} exp(V_{p_t}) + sum_{m in B} exp(V_m)),
//   one stage per position; a group with nothing below it contributes 1. A
//   ranking without ties is the product over positions of the multinomial-logit
//   choice of the ranked alternative among those not yet ranked (Beggs, Cardell
//   and Hausman 1981).
//   Stability across stages (split > 0, rank_test()): a stage at position
//   t > split uses a second parameter block, theta = (theta_early, theta_late).
//
// Parameter vector (one block, twice when split > 0):
//   beta  (p * (J-1))  alternative-major, as in mnl_ll.cpp;   gamma (q)
//
// Derivatives are formed in V (one J-vector per block) and mapped to theta
// through the design rows. A stage choosing c from S with probabilities P
// contributes e_c - P to the gradient and -(diag(P) - P P') to the Hessian; the
// orders p of a tie group, with log-probabilities l_p, V-gradients g_p and
// V-Hessians H_p, weights omega_p = exp(l_p) / sum_q exp(l_q), contribute
// sum omega_p g_p and sum omega_p (H_p + g_p g_p') - gbar gbar' (the
// derivatives of the log of a sum). The orders are accumulated with a running
// log-sum-exp. Validated against numDeriv and against mnl() on the exploded
// data in tests/testthat/test-rologit.R.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
using namespace Rcpp;

namespace {

struct RSpec {
  int n, J, p, q, base, K0, nb, split, K, M;   // M = nb * J, the size of the V-space
  std::vector<int> bpos;
  int gpos;
};

RSpec make_rspec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
                 const NumericMatrix& off, const IntegerMatrix& avail, const IntegerMatrix& R,
                 int J, int base, int split) {
  RSpec s;
  s.n = R.nrow(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.base = base; s.split = split;
  s.bpos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  s.K0 = pos; s.nb = (split > 0) ? 2 : 1; s.K = s.nb * s.K0; s.M = s.nb * J;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (X.nrow() != s.n) stop("internal: design rows do not match the rankings.");
  if (W.nrow() != s.n * J) stop("internal: W must have n * J rows.");
  if (R.ncol() != J || off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J)
    stop("internal: rank, offset, or availability dimensions do not match.");
  return s;
}

struct Work {
  std::vector<double> V, P, g, H, gp, Hp, G, HH;
  std::vector<int> cand, below, T;
  std::vector< std::pair<int, int> > ranked;
  std::vector<int> unranked;
  void init(const RSpec& s) {
    V.assign(s.M, 0.0); P.assign(s.J, 0.0); g.assign(s.M, 0.0); H.assign(s.M * s.M, 0.0);
    gp.assign(s.M, 0.0); Hp.assign(s.M * s.M, 0.0); G.assign(s.M, 0.0); HH.assign(s.M * s.M, 0.0);
  }
};

void utilities(int i, const RSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
               const NumericMatrix& off, const IntegerMatrix& avail, Work& w) {
  const int J = s.J;
  for (int b = 0; b < s.nb; ++b) {
    const int o = b * s.K0;
    for (int j = 0; j < J; ++j) {
      if (!avail(i, j)) { w.V[b * J + j] = R_NegInf; continue; }
      double v = off(i, j);
      if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) v += X(i, c) * theta[o + s.bpos[j] + c];
      for (int c = 0; c < s.q; ++c) v += W(i * J + j, c) * theta[o + s.gpos + c];
      w.V[b * J + j] = v;
    }
  }
}

// one stage: alternative c chosen among w.cand with block b; returns the
// log-probability and adds the V-derivatives to g (and H when non-null)
double stage(int c, int b, const RSpec& s, Work& w, double* g, double* H) {
  const int J = s.J, M = s.M;
  const double* V = &w.V[b * J];
  double vmax = R_NegInf;
  for (size_t a = 0; a < w.cand.size(); ++a) vmax = std::max(vmax, V[w.cand[a]]);
  double den = 0.0;
  for (size_t a = 0; a < w.cand.size(); ++a) den += std::exp(V[w.cand[a]] - vmax);
  const double lden = vmax + std::log(den);
  for (size_t a = 0; a < w.cand.size(); ++a) w.P[w.cand[a]] = std::exp(V[w.cand[a]] - lden);
  if (g) for (size_t a = 0; a < w.cand.size(); ++a) { const int j = w.cand[a]; g[b * J + j] += (j == c ? 1.0 : 0.0) - w.P[j]; }
  if (H) for (size_t a = 0; a < w.cand.size(); ++a) {
    const int j = w.cand[a];
    for (size_t e = 0; e < w.cand.size(); ++e) {
      const int m = w.cand[e];
      H[(b * J + j) * M + b * J + m] -= w.P[j] * ((j == m ? 1.0 : 0.0) - w.P[m]);
    }
  }
  return V[c] - lden;
}

inline int block_at(const RSpec& s, int position) { return (s.split > 0 && position > s.split) ? 1 : 0; }

// the log-likelihood of case i; with want_g / want_h the V-gradient and
// V-Hessian of the case are left in w.g / w.H
double eval_case(int i, const RSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
                 const NumericMatrix& off, const IntegerMatrix& avail, const IntegerMatrix& R, Work& w,
                 bool want_g, bool want_h) {
  const int J = s.J, M = s.M;
  utilities(i, s, theta, X, W, off, avail, w);
  if (want_g) std::fill(w.g.begin(), w.g.end(), 0.0);
  if (want_h) std::fill(w.H.begin(), w.H.end(), 0.0);
  w.ranked.clear(); w.unranked.clear();
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    if (R(i, j) > 0) w.ranked.push_back(std::make_pair(R(i, j), j)); else w.unranked.push_back(j);
  }
  std::sort(w.ranked.begin(), w.ranked.end());
  double ll = 0.0;
  int pos = 0;                                     // alternatives placed above the current group
  size_t a = 0;
  while (a < w.ranked.size()) {
    size_t e = a;
    while (e < w.ranked.size() && w.ranked[e].first == w.ranked[a].first) ++e;
    const int d = (int) (e - a);
    w.below.clear();
    for (size_t t = e; t < w.ranked.size(); ++t) w.below.push_back(w.ranked[t].second);
    for (size_t t = 0; t < w.unranked.size(); ++t) w.below.push_back(w.unranked[t]);
    if (w.below.empty()) break;                    // nothing below: the group's term is 1
    w.T.clear();
    for (size_t t = a; t < e; ++t) w.T.push_back(w.ranked[t].second);
    if (d == 1) {
      w.cand = w.below; w.cand.push_back(w.T[0]);
      ll += stage(w.T[0], block_at(s, pos + 1), s, w, want_g ? w.g.data() : NULL, want_h ? w.H.data() : NULL);
    } else {
      if (d > 7) stop("internal: a tie group of %d alternatives (at most 7 are enumerated).", d);
      std::sort(w.T.begin(), w.T.end());
      // the Hessian of the log of a sum needs every order's gradient too
      const bool need_g = want_g || want_h;
      double mx = R_NegInf, S = 0.0;
      if (need_g) std::fill(w.G.begin(), w.G.end(), 0.0);
      if (want_h) std::fill(w.HH.begin(), w.HH.end(), 0.0);
      do {
        if (need_g) std::fill(w.gp.begin(), w.gp.end(), 0.0);
        if (want_h) std::fill(w.Hp.begin(), w.Hp.end(), 0.0);
        double lp = 0.0;
        for (int st = 0; st < d; ++st) {
          w.cand = w.below;
          for (int t = st; t < d; ++t) w.cand.push_back(w.T[t]);
          lp += stage(w.T[st], block_at(s, pos + st + 1), s, w, need_g ? w.gp.data() : NULL, want_h ? w.Hp.data() : NULL);
        }
        if (lp > mx) {                             // rescale the running sums to the new maximum
          const double sc = (mx == R_NegInf) ? 0.0 : std::exp(mx - lp);
          S *= sc;
          if (need_g) for (int u = 0; u < M; ++u) w.G[u] *= sc;
          if (want_h) for (int u = 0; u < M * M; ++u) w.HH[u] *= sc;
          mx = lp;
        }
        const double om = std::exp(lp - mx);
        S += om;
        if (need_g) for (int u = 0; u < M; ++u) w.G[u] += om * w.gp[u];
        if (want_h) for (int u = 0; u < M; ++u) for (int v = 0; v < M; ++v)
          w.HH[u * M + v] += om * (w.Hp[u * M + v] + w.gp[u] * w.gp[v]);
      } while (std::next_permutation(w.T.begin(), w.T.end()));
      ll += mx + std::log(S);
      if (need_g) for (int u = 0; u < M; ++u) w.G[u] /= S;
      if (want_g) for (int u = 0; u < M; ++u) w.g[u] += w.G[u];
      if (want_h) for (int u = 0; u < M; ++u) for (int v = 0; v < M; ++v)
        w.H[u * M + v] += w.HH[u * M + v] / S - w.G[u] * w.G[v];
    }
    pos += d;
    a = e;
  }
  return ll;
}

// the theta-score of case i from its V-gradient, into out (overwritten)
void score_theta(int i, const RSpec& s, const NumericMatrix& X, const NumericMatrix& W, const IntegerMatrix& avail,
                 const Work& w, std::vector<double>& out) {
  std::fill(out.begin(), out.end(), 0.0);
  const int J = s.J;
  for (int b = 0; b < s.nb; ++b) {
    const int o = b * s.K0;
    for (int j = 0; j < J; ++j) {
      if (!avail(i, j)) continue;
      const double r = w.g[b * J + j];
      if (r == 0.0) continue;
      if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) out[o + s.bpos[j] + c] += r * X(i, c);
      for (int c = 0; c < s.q; ++c) out[o + s.gpos + c] += r * W(i * J + j, c);
    }
  }
}

} // namespace

// Per-case log-likelihood contributions (unweighted).
// [[Rcpp::export]]
NumericVector pc_ro_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off,
                            IntegerMatrix avail, IntegerMatrix R, int J, int base, int split) {
  RSpec s = make_rspec(theta, X, W, off, avail, R, J, base, split);
  Work w; w.init(s);
  NumericVector out(s.n);
  for (int i = 0; i < s.n; ++i) out[i] = eval_case(i, s, theta, X, W, off, avail, R, w, false, false);
  return out;
}

// Gradient of the weighted log-likelihood.
// [[Rcpp::export]]
NumericVector pc_ro_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off,
                             IntegerMatrix avail, IntegerMatrix R, NumericVector wt, int J, int base, int split) {
  RSpec s = make_rspec(theta, X, W, off, avail, R, J, base, split);
  if (wt.size() != s.n) stop("internal: weights do not match the cases.");
  Work w; w.init(s);
  std::vector<double> sc(s.K);
  NumericVector g(s.K);
  for (int i = 0; i < s.n; ++i) {
    if (wt[i] == 0.0) continue;
    eval_case(i, s, theta, X, W, off, avail, R, w, true, false);
    score_theta(i, s, X, W, avail, w, sc);
    for (int c = 0; c < s.K; ++c) g[c] += wt[i] * sc[c];
  }
  return g;
}

// Per-case scores (unweighted), n x K.
// [[Rcpp::export]]
NumericMatrix pc_ro_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off,
                               IntegerMatrix avail, IntegerMatrix R, int J, int base, int split) {
  RSpec s = make_rspec(theta, X, W, off, avail, R, J, base, split);
  Work w; w.init(s);
  std::vector<double> sc(s.K);
  NumericMatrix S(s.n, s.K);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, off, avail, R, w, true, false);
    score_theta(i, s, X, W, avail, w, sc);
    for (int c = 0; c < s.K; ++c) S(i, c) = sc[c];
  }
  return S;
}

// Hessian of the weighted log-likelihood (K x K): the V-Hessian of each case
// mapped through its sparse design rows D_(b,j) = dV_j^(b) / dtheta.
// [[Rcpp::export]]
NumericMatrix pc_ro_hess_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off,
                             IntegerMatrix avail, IntegerMatrix R, NumericVector wt, int J, int base, int split) {
  RSpec s = make_rspec(theta, X, W, off, avail, R, J, base, split);
  if (wt.size() != s.n) stop("internal: weights do not match the cases.");
  Work w; w.init(s);
  const int M = s.M;
  std::vector< std::vector< std::pair<int, double> > > D(M);
  NumericMatrix Hm(s.K, s.K);
  for (int i = 0; i < s.n; ++i) {
    const double wi = wt[i];
    if (wi == 0.0) continue;
    eval_case(i, s, theta, X, W, off, avail, R, w, false, true);
    for (int b = 0; b < s.nb; ++b) {
      const int o = b * s.K0;
      for (int j = 0; j < J; ++j) {
        std::vector< std::pair<int, double> >& Dj = D[b * J + j];
        Dj.clear();
        if (!avail(i, j)) continue;
        if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) { const double v = X(i, c); if (v != 0.0) Dj.push_back(std::make_pair(o + s.bpos[j] + c, v)); }
        for (int c = 0; c < s.q; ++c) { const double v = W(i * J + j, c); if (v != 0.0) Dj.push_back(std::make_pair(o + s.gpos + c, v)); }
      }
    }
    for (int u = 0; u < M; ++u) {
      if (D[u].empty()) continue;
      for (int v = 0; v < M; ++v) {
        if (D[v].empty()) continue;
        const double h = wi * w.H[u * M + v];
        if (h == 0.0) continue;
        for (size_t a = 0; a < D[u].size(); ++a)
          for (size_t b = 0; b < D[v].size(); ++b)
            Hm(D[u][a].first, D[v][b].first) += h * D[u][a].second * D[v][b].second;
      }
    }
  }
  return Hm;
}

// First-choice probabilities (the multinomial logit over the choice set, with
// the first parameter block) and the utilities V (n x J, -Inf where unavailable).
// [[Rcpp::export]]
List pc_ro_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off,
                     IntegerMatrix avail, int J, int base) {
  const int n = X.nrow(), p = X.ncol(), q = W.ncol();
  std::vector<int> bpos(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { bpos[j] = pos; pos += p; }
  const int gpos = pos;
  if (theta.size() < gpos + q) stop("internal: theta is too short.");
  NumericMatrix P(n, J), Vm(n, J);
  std::vector<double> V(J);
  for (int i = 0; i < n; ++i) {
    double vmax = R_NegInf;
    for (int j = 0; j < J; ++j) {
      if (!avail(i, j)) { V[j] = R_NegInf; continue; }
      double v = off(i, j);
      if (bpos[j] >= 0) for (int c = 0; c < p; ++c) v += X(i, c) * theta[bpos[j] + c];
      for (int c = 0; c < q; ++c) v += W(i * J + j, c) * theta[gpos + c];
      V[j] = v; vmax = std::max(vmax, v);
    }
    double den = 0.0;
    for (int j = 0; j < J; ++j) if (V[j] != R_NegInf) den += std::exp(V[j] - vmax);
    for (int j = 0; j < J; ++j) { P(i, j) = (V[j] == R_NegInf) ? 0.0 : std::exp(V[j] - vmax) / den; Vm(i, j) = V[j]; }
  }
  return List::create(_["prob"] = P, _["V"] = Vm);
}
