// nlogit_ll.cpp -- the two-level nested logit (McFadden 1978) in the
// utility-maximization-consistent ("RU2") form, the one mlogit and Stata's
// nlogit (RUM-consistent) use:
//
//   V_ij as in mnl_ll.cpp (case-specific x with alternative-specific beta,
//   alternative-specific w with generic gamma, offsets, availability)
//   nest m with inclusive-value parameter lambda_m:
//     I_im   = log sum_{k in B_m, available} exp(V_ik / lambda_m)
//     P(j | m) = exp(V_ij / lambda_m - I_im)
//     P(m)     = exp(lambda_m I_im) / sum_l exp(lambda_l I_il)     (over nests with an available member)
//     P(j)     = P(j | m(j)) P(m(j))
//
// Parameter vector theta: beta (p (J-1)), gamma (q), then ell, the log of the
// free inclusive-value parameters: lmap[m] is the index into the ell block for
// nest m (-1 when lambda_m is fixed at 1: a singleton nest, whose lambda does
// not enter the likelihood). With equal lambdas every nest maps to the same
// entry.
//
// Per-case derivatives (y the chosen alternative, in nest m_y):
//   d l / d V_k   = [1{k = y} - P(k|m_y)] / lambda_{m_y} + P(k|m_y) (1 - P(m_y))      k in B_{m_y}
//                 = - P(l) P(k|l)                                                     k in B_l, l != m_y
//   d l / d lambda_{m_y} = (Ubar_{m_y} - V_y) / lambda^2 + (1 - P(m_y)) (I_{m_y} - Ubar_{m_y} / lambda)
//   d l / d lambda_l     = - P(l) (I_l - Ubar_l / lambda_l)                              l != m_y
//   with Ubar_m = sum_{k in B_m} P(k|m) V_k;  d l / d ell = lambda d l / d lambda.
// Validated against numDeriv in tests/testthat/test-families.R. The Hessian is
// taken by finite differences of this gradient (R side).
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

namespace {

struct NLSpec {
  int n, J, p, q, nnest, nell, base, K;
  std::vector<int> bpos;
  int gpos, lpos;
};

NLSpec nl_spec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W, const IntegerVector& y,
               int J, int base, const IntegerVector& nest, const IntegerVector& lmap) {
  NLSpec s;
  s.n = y.size(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.base = base;
  s.nnest = lmap.size();
  int nell = 0; for (int m = 0; m < s.nnest; ++m) if (lmap[m] + 1 > nell) nell = lmap[m] + 1;
  s.nell = nell;
  s.bpos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  s.lpos = pos; pos += s.nell;
  s.K = pos;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (nest.size() != J) stop("internal: nest must have one entry per alternative.");
  if (X.nrow() != s.n || W.nrow() != s.n * J) stop("internal: design rows do not match the response.");
  return s;
}

inline double lse(const std::vector<double>& v, const std::vector<int>& use) {
  double m = R_NegInf;
  for (size_t i = 0; i < use.size(); ++i) if (v[use[i]] > m) m = v[use[i]];
  if (m == R_NegInf) return R_NegInf;
  double s = 0.0;
  for (size_t i = 0; i < use.size(); ++i) s += std::exp(v[use[i]] - m);
  return m + std::log(s);
}

// one case: fills V, the conditional probabilities Pc (P(k|nest)), nest
// probabilities Pn, inclusive values I, Ubar, and the log-likelihood
struct NLCase {
  std::vector<double> V, Pc, Pn, I, Ubar, lam, P;
  std::vector<int> nonempty;
  double lli;
};

void nl_case(int i, const NLSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
             const NumericMatrix& off, const IntegerMatrix& avail, const IntegerVector& nest, const IntegerVector& lmap,
             int yi, NLCase& c) {
  const int J = s.J, M = s.nnest;
  c.V.assign(J, R_NegInf); c.Pc.assign(J, 0.0); c.P.assign(J, 0.0);
  c.Pn.assign(M, 0.0); c.I.assign(M, R_NegInf); c.Ubar.assign(M, 0.0); c.lam.assign(M, 1.0);
  for (int m = 0; m < M; ++m) if (lmap[m] >= 0) c.lam[m] = std::exp(theta[s.lpos + lmap[m]]);
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    double v = off(i, j);
    if (s.bpos[j] >= 0) for (int cc = 0; cc < s.p; ++cc) v += X(i, cc) * theta[s.bpos[j] + cc];
    for (int cc = 0; cc < s.q; ++cc) v += W(i * J + j, cc) * theta[s.gpos + cc];
    c.V[j] = v;
  }
  // scaled utilities per nest
  std::vector<double> a(J, R_NegInf);
  c.nonempty.clear();
  for (int m = 0; m < M; ++m) {
    std::vector<int> mem;
    for (int j = 0; j < J; ++j) if (nest[j] == m && avail(i, j)) { mem.push_back(j); a[j] = c.V[j] / c.lam[m]; }
    if (mem.empty()) continue;
    c.nonempty.push_back(m);
    c.I[m] = lse(a, mem);
    double ub = 0.0;
    for (size_t t = 0; t < mem.size(); ++t) { int j = mem[t]; c.Pc[j] = std::exp(a[j] - c.I[m]); ub += c.Pc[j] * c.V[j]; }
    c.Ubar[m] = ub;
  }
  // nest probabilities. B_m = lambda_m I_m is split as Vmax_m + lambda_m log n_m + C_m
  // with C_m = lambda_m log1p(mean_k expm1(-(Vmax_m - V_k) / lambda_m)) = O(1), and the
  // nests are compared part by part: summing first would round the O(1) part
  // away once lambda is large (with equal lambdas and nest sizes the O(lambda)
  // parts cancel exactly and the O(1) parts decide; the summed form returned
  // log P = 0 for every case at lambda = e^40, a false optimum)
  if (c.nonempty.empty()) Rcpp::stop("internal: a case has no available alternative.");
  std::vector<double> Bv(M, 0.0), Bl(M, 0.0), d(M, R_NegInf);
  for (size_t t = 0; t < c.nonempty.size(); ++t) {
    const int m = c.nonempty[t];
    double vmax = R_NegInf; int nm = 0;
    for (int j = 0; j < J; ++j) if (nest[j] == m && avail(i, j)) { ++nm; if (c.V[j] > vmax) vmax = c.V[j]; }
    double sm = 0.0;
    for (int j = 0; j < J; ++j) if (nest[j] == m && avail(i, j)) sm += std::expm1(-(vmax - c.V[j]) / c.lam[m]);
    Bv[m] = vmax + c.lam[m] * std::log1p(sm / nm);
    Bl[m] = c.lam[m] * std::log((double) nm);
  }
  const int ref = c.nonempty[0];
  for (size_t t = 0; t < c.nonempty.size(); ++t) { int m = c.nonempty[t]; d[m] = (Bv[m] - Bv[ref]) + (Bl[m] - Bl[ref]); }
  const double lseD = lse(d, c.nonempty);
  for (size_t t = 0; t < c.nonempty.size(); ++t) { int m = c.nonempty[t]; c.Pn[m] = std::exp(d[m] - lseD); }
  for (int j = 0; j < J; ++j) if (avail(i, j)) c.P[j] = c.Pc[j] * c.Pn[nest[j]];
  const int my = nest[yi];
  c.lli = (a[yi] - c.I[my]) + (d[my] - lseD);
}

void nl_score(int i, const NLSpec& s, const NLCase& c, const NumericMatrix& X, const NumericMatrix& W,
              const IntegerMatrix& avail, const IntegerVector& nest, const IntegerVector& lmap, int yi, std::vector<double>& out) {
  std::fill(out.begin(), out.end(), 0.0);
  const int J = s.J, my = nest[yi];
  for (int k = 0; k < J; ++k) {
    if (!avail(i, k)) continue;
    const int l = nest[k];
    double dV;
    if (l == my) dV = ((k == yi ? 1.0 : 0.0) - c.Pc[k]) / c.lam[my] + c.Pc[k] * (1.0 - c.Pn[my]);
    else dV = -c.Pn[l] * c.Pc[k];
    if (s.bpos[k] >= 0) for (int cc = 0; cc < s.p; ++cc) out[s.bpos[k] + cc] += dV * X(i, cc);
    for (int cc = 0; cc < s.q; ++cc) out[s.gpos + cc] += dV * W(i * J + k, cc);
  }
  for (size_t t = 0; t < c.nonempty.size(); ++t) {
    const int m = c.nonempty[t];
    if (lmap[m] < 0) continue;
    const double lam = c.lam[m];
    double dl;
    if (m == my) dl = (c.Ubar[m] - c.V[yi]) / (lam * lam) + (1.0 - c.Pn[m]) * (c.I[m] - c.Ubar[m] / lam);
    else dl = -c.Pn[m] * (c.I[m] - c.Ubar[m] / lam);
    out[s.lpos + lmap[m]] += dl * lam;                 // d/d ell = lambda d/d lambda
  }
}

} // namespace

// [[Rcpp::export]]
NumericVector pc_nl_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                            IntegerVector y, int J, int base, IntegerVector nest, IntegerVector lmap) {
  NLSpec s = nl_spec(theta, X, W, y, J, base, nest, lmap);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  NLCase c; NumericVector out(s.n);
  for (int i = 0; i < s.n; ++i) { nl_case(i, s, theta, X, W, off, avail, nest, lmap, y[i], c); out[i] = c.lli; }
  return out;
}

// [[Rcpp::export]]
NumericVector pc_nl_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                             IntegerVector y, NumericVector w, int J, int base, IntegerVector nest, IntegerVector lmap) {
  NLSpec s = nl_spec(theta, X, W, y, J, base, nest, lmap);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  if (w.size() != s.n) stop("internal: weights do not match the cases.");
  NLCase c; std::vector<double> sc(s.K); NumericVector g(s.K);
  for (int i = 0; i < s.n; ++i) {
    nl_case(i, s, theta, X, W, off, avail, nest, lmap, y[i], c);
    nl_score(i, s, c, X, W, avail, nest, lmap, y[i], sc);
    for (int k = 0; k < s.K; ++k) g[k] += w[i] * sc[k];
  }
  return g;
}

// [[Rcpp::export]]
NumericMatrix pc_nl_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                               IntegerVector y, int J, int base, IntegerVector nest, IntegerVector lmap) {
  NLSpec s = nl_spec(theta, X, W, y, J, base, nest, lmap);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  NLCase c; std::vector<double> sc(s.K); NumericMatrix S(s.n, s.K);
  for (int i = 0; i < s.n; ++i) {
    nl_case(i, s, theta, X, W, off, avail, nest, lmap, y[i], c);
    nl_score(i, s, c, X, W, avail, nest, lmap, y[i], sc);
    for (int k = 0; k < s.K; ++k) S(i, k) = sc[k];
  }
  return S;
}

// probabilities P(j) (n x J), the conditional P(j | nest) and the utilities V
// [[Rcpp::export]]
List pc_nl_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                     IntegerVector y, int J, int base, IntegerVector nest, IntegerVector lmap) {
  NLSpec s = nl_spec(theta, X, W, y, J, base, nest, lmap);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  NLCase c; NumericMatrix P(s.n, J), Pc(s.n, J), Vm(s.n, J);
  for (int i = 0; i < s.n; ++i) {
    int yi = y[i];
    if (yi < 0 || yi >= J || !avail(i, yi)) { yi = 0; while (yi < J && !avail(i, yi)) ++yi; if (yi == J) yi = 0; }
    nl_case(i, s, theta, X, W, off, avail, nest, lmap, yi, c);
    for (int j = 0; j < J; ++j) { P(i, j) = c.P[j]; Pc(i, j) = c.Pc[j]; Vm(i, j) = c.V[j]; }
  }
  return List::create(_["prob"] = P, _["prob_conditional"] = Pc, _["V"] = Vm);
}
