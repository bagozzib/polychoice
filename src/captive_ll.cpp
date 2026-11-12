// captive_ll.cpp -- the likelihood engine of imnl() with several inflated
// (captive) alternatives: the multinomial logit with covariate-dependent
// captive shares, i.e. the dogit of Gaudry and Dagenais (1979) with captivity
// parameters theta_ic = exp(-z_i' delta_c). imnl() with one inflated
// alternative keeps mnl_ll.cpp; at |C| = 1 this kernel is the same model and
// is checked against that one in tests/testthat/test-captive.R.
//
// Model (case i = 0..n-1, alternatives j = 0..J-1, base b, availability a_ij;
// C the captive alternatives):
//   V_ij, P_ij as in mnl_ll.cpp (the choice regime)
//   eta_ic = z_i' delta_c for c in C with a_ic = 1 (a captive alternative
//   outside the case's choice set is not one of its regimes: the dogit sums
//   the captivity terms over the choice set)
//   pi_i0 = 1 / D_i,  pi_ic = exp(-eta_ic) / D_i,  D_i = 1 + sum_c exp(-eta_ic)
//   P*_ij = pi_i0 P_ij + 1{j in C} pi_ij,   l_i = log P*_{i, y_i}
// so eta_ic = log(pi_i0 / pi_ic): a positive coefficient makes captivity to c
// less likely relative to choosing, the reading of imnl()'s split equation.
//
// Parameter vector theta: beta (p * (J-1), alternative-major as in
// mnl_ll.cpp), gamma (q), delta (r * |C|: the r split coefficients of each
// captive alternative in the order of C).
//
// Derivatives (DESIGN.md, Phase 3 A). With d_j = 1{j = y} - P_j, the posterior
// of the choice regime tau0 = pi_0 P_y / P*_y and of captivity
// rho_c = 1{y = c} pi_c / P*_y:
//   dl/dV_j            = tau0 d_j
//   dl/deta_c          = pi_c - rho_c
//   d2l/dV_j dV_m      = tau0 (1 - tau0) d_j d_m - tau0 P_j (1{j = m} - P_m)
//   d2l/dV_j deta_c    = tau0 rho_c d_j
//   d2l/deta_c deta_e  = pi_c (pi_e - 1{c = e}) - rho_c (rho_e - 1{c = e})
// The regime probabilities are formed on the log scale (log-sum-exp over the
// case's regimes) and log P*_y as a two-term log-sum-exp.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
using namespace Rcpp;

namespace {

struct CSpec {
  int n, J, p, q, r, base, nc, K;
  std::vector<int> cap;      // the captive alternatives (0-based), in the order of their delta blocks
  std::vector<int> capidx;   // capidx[j]: position of j in cap, -1 if j is not captive
  std::vector<int> bpos;
  int gpos, dpos;
};

CSpec make_cspec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W, const NumericMatrix& Z,
                 const NumericMatrix& off, const IntegerMatrix& avail, const IntegerVector& y,
                 int J, int base, const IntegerVector& cap) {
  CSpec s;
  s.n = y.size(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.r = Z.ncol(); s.base = base; s.nc = cap.size();
  if (s.nc < 1) stop("internal: no captive alternative.");
  s.cap.assign(cap.begin(), cap.end());
  s.capidx.assign(J, -1);
  for (int t = 0; t < s.nc; ++t) {
    if (s.cap[t] < 0 || s.cap[t] >= J || s.capidx[s.cap[t]] >= 0) stop("internal: invalid captive alternatives.");
    s.capidx[s.cap[t]] = t;
  }
  s.bpos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  s.dpos = pos; pos += s.r * s.nc;
  s.K = pos;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (X.nrow() != s.n || Z.nrow() != s.n) stop("internal: design rows do not match the response.");
  if (W.nrow() != s.n * J) stop("internal: W must have n * J rows.");
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  return s;
}

inline double lse2(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  const double m = std::max(a, b);
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

struct CEval {
  std::vector<double> V, P, Pstar, pi, rho, d;   // pi, rho: one per captive (0 where unavailable)
  double pi0, logpi0, tau0, lli;
  void init(const CSpec& s) {
    V.assign(s.J, 0.0); P.assign(s.J, 0.0); Pstar.assign(s.J, 0.0); d.assign(s.J, 0.0);
    pi.assign(s.nc, 0.0); rho.assign(s.nc, 0.0);
  }
};

void eval_case(int i, const CSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
               const NumericMatrix& Z, const NumericMatrix& off, const IntegerMatrix& avail, int yi, CEval& e) {
  const int J = s.J;
  double vmax = R_NegInf;
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) { e.V[j] = R_NegInf; continue; }
    double v = off(i, j);
    if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) v += X(i, c) * theta[s.bpos[j] + c];
    for (int c = 0; c < s.q; ++c) v += W(i * J + j, c) * theta[s.gpos + c];
    e.V[j] = v;
    if (v > vmax) vmax = v;
  }
  double den = 0.0;
  for (int j = 0; j < J; ++j) if (e.V[j] != R_NegInf) den += std::exp(e.V[j] - vmax);
  const double logden = vmax + std::log(den);
  for (int j = 0; j < J; ++j) e.P[j] = (e.V[j] == R_NegInf) ? 0.0 : std::exp(e.V[j] - logden);
  const double logPy = e.V[yi] - logden;
  // the regimes: log D = log(1 + sum_c exp(-eta_c)) over the available captives
  std::vector<double> a(s.nc, R_NegInf);
  double logD = 0.0;
  for (int t = 0; t < s.nc; ++t) {
    if (!avail(i, s.cap[t])) continue;
    double eta = 0.0;
    for (int c = 0; c < s.r; ++c) eta += Z(i, c) * theta[s.dpos + t * s.r + c];
    a[t] = -eta;
    logD = lse2(logD, a[t]);
  }
  e.logpi0 = -logD; e.pi0 = std::exp(e.logpi0);
  for (int t = 0; t < s.nc; ++t) e.pi[t] = (a[t] == R_NegInf) ? 0.0 : std::exp(a[t] - logD);
  for (int j = 0; j < J; ++j) {
    const int t = s.capidx[j];
    e.Pstar[j] = e.pi0 * e.P[j] + ((t >= 0) ? e.pi[t] : 0.0);
  }
  const int ty = s.capidx[yi];
  double logPstar_y = e.logpi0 + logPy;
  if (ty >= 0 && a[ty] != R_NegInf) logPstar_y = lse2(logPstar_y, a[ty] - logD);
  e.lli = logPstar_y;
  e.tau0 = std::exp(e.logpi0 + logPy - logPstar_y);
  for (int t = 0; t < s.nc; ++t) e.rho[t] = (t == ty && a[t] != R_NegInf) ? std::exp(a[t] - logD - logPstar_y) : 0.0;
  for (int j = 0; j < J; ++j) e.d[j] = (j == yi ? 1.0 : 0.0) - e.P[j];
}

void score_case(int i, const CSpec& s, const CEval& e, const NumericMatrix& X, const NumericMatrix& W,
                const NumericMatrix& Z, const IntegerMatrix& avail, std::vector<double>& out) {
  std::fill(out.begin(), out.end(), 0.0);
  const int J = s.J;
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    const double rj = e.tau0 * e.d[j];
    if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) out[s.bpos[j] + c] += rj * X(i, c);
    for (int c = 0; c < s.q; ++c) out[s.gpos + c] += rj * W(i * J + j, c);
  }
  for (int t = 0; t < s.nc; ++t) {
    const double u = e.pi[t] - e.rho[t];
    if (u == 0.0) continue;
    for (int c = 0; c < s.r; ++c) out[s.dpos + t * s.r + c] = u * Z(i, c);
  }
}

} // namespace

// Per-case log-likelihood contributions log P*_{i, y_i} (unweighted).
// [[Rcpp::export]]
NumericVector pc_cap_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                             NumericMatrix off, IntegerMatrix avail, IntegerVector y, int J, int base, IntegerVector cap) {
  CSpec s = make_cspec(theta, X, W, Z, off, avail, y, J, base, cap);
  CEval e; e.init(s);
  NumericVector out(s.n);
  for (int i = 0; i < s.n; ++i) { eval_case(i, s, theta, X, W, Z, off, avail, y[i], e); out[i] = e.lli; }
  return out;
}

// Gradient of the weighted log-likelihood.
// [[Rcpp::export]]
NumericVector pc_cap_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                              NumericMatrix off, IntegerMatrix avail, IntegerVector y, NumericVector w,
                              int J, int base, IntegerVector cap) {
  CSpec s = make_cspec(theta, X, W, Z, off, avail, y, J, base, cap);
  if (w.size() != s.n) stop("internal: weights do not match the cases.");
  CEval e; e.init(s);
  std::vector<double> sc(s.K);
  NumericVector g(s.K);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, Z, off, avail, y[i], e);
    score_case(i, s, e, X, W, Z, avail, sc);
    for (int c = 0; c < s.K; ++c) g[c] += w[i] * sc[c];
  }
  return g;
}

// Per-case scores (unweighted), n x K.
// [[Rcpp::export]]
NumericMatrix pc_cap_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                                NumericMatrix off, IntegerMatrix avail, IntegerVector y, int J, int base, IntegerVector cap) {
  CSpec s = make_cspec(theta, X, W, Z, off, avail, y, J, base, cap);
  CEval e; e.init(s);
  std::vector<double> sc(s.K);
  NumericMatrix S(s.n, s.K);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, Z, off, avail, y[i], e);
    score_case(i, s, e, X, W, Z, avail, sc);
    for (int c = 0; c < s.K; ++c) S(i, c) = sc[c];
  }
  return S;
}

// Hessian of the weighted log-likelihood (K x K).
// [[Rcpp::export]]
NumericMatrix pc_cap_hess_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                              NumericMatrix off, IntegerMatrix avail, IntegerVector y, NumericVector w,
                              int J, int base, IntegerVector cap) {
  CSpec s = make_cspec(theta, X, W, Z, off, avail, y, J, base, cap);
  if (w.size() != s.n) stop("internal: weights do not match the cases.");
  CEval e; e.init(s);
  std::vector< std::vector< std::pair<int, double> > > D(J);
  NumericMatrix H(s.K, s.K);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, Z, off, avail, y[i], e);
    const double wi = w[i];
    if (wi == 0.0) continue;
    for (int j = 0; j < J; ++j) {
      D[j].clear();
      if (!avail(i, j)) continue;
      if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) { const double v = X(i, c); if (v != 0.0) D[j].push_back(std::make_pair(s.bpos[j] + c, v)); }
      for (int c = 0; c < s.q; ++c) { const double v = W(i * J + j, c); if (v != 0.0) D[j].push_back(std::make_pair(s.gpos + c, v)); }
    }
    const double t0 = e.tau0, a1 = t0 * (1.0 - t0);
    // V-V block
    for (int j = 0; j < J; ++j) {
      if (D[j].empty()) continue;
      for (int m = 0; m < J; ++m) {
        if (D[m].empty()) continue;
        const double h = wi * (a1 * e.d[j] * e.d[m] - t0 * e.P[j] * ((j == m ? 1.0 : 0.0) - e.P[m]));
        if (h == 0.0) continue;
        for (size_t u = 0; u < D[j].size(); ++u)
          for (size_t v = 0; v < D[m].size(); ++v)
            H(D[j][u].first, D[m][v].first) += h * D[j][u].second * D[m][v].second;
      }
    }
    // V-eta cross terms: tau0 rho_t d_j (only the chosen captive has rho != 0)
    for (int t = 0; t < s.nc; ++t) {
      if (e.rho[t] == 0.0) continue;
      for (int j = 0; j < J; ++j) {
        if (D[j].empty()) continue;
        const double hv = wi * t0 * e.rho[t] * e.d[j];
        if (hv == 0.0) continue;
        for (size_t u = 0; u < D[j].size(); ++u)
          for (int c = 0; c < s.r; ++c) {
            const double v = hv * D[j][u].second * Z(i, c);
            H(D[j][u].first, s.dpos + t * s.r + c) += v;
            H(s.dpos + t * s.r + c, D[j][u].first) += v;
          }
      }
    }
    // eta-eta blocks
    for (int t = 0; t < s.nc; ++t) {
      if (e.pi[t] == 0.0 && e.rho[t] == 0.0) continue;
      for (int t2 = 0; t2 < s.nc; ++t2) {
        const double kd = (t == t2) ? 1.0 : 0.0;
        const double h = wi * (e.pi[t] * (e.pi[t2] - kd) - e.rho[t] * (e.rho[t2] - kd));
        if (h == 0.0) continue;
        for (int c = 0; c < s.r; ++c)
          for (int c2 = 0; c2 < s.r; ++c2)
            H(s.dpos + t * s.r + c, s.dpos + t2 * s.r + c2) += h * Z(i, c) * Z(i, c2);
      }
    }
  }
  return H;
}

// Probabilities: P* (observed categories), P (choice regime), the regime
// probabilities (n x (1 + |C|): the choice regime, then captivity to each
// alternative of C, 0 where it is unavailable), phi = P(choice regime), and V.
// `y` is only used for its length; pass any valid alternative.
// [[Rcpp::export]]
List pc_cap_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                      NumericMatrix off, IntegerMatrix avail, IntegerVector y, int J, int base, IntegerVector cap) {
  CSpec s = make_cspec(theta, X, W, Z, off, avail, y, J, base, cap);
  CEval e; e.init(s);
  NumericMatrix Pstar(s.n, J), P(s.n, J), Vm(s.n, J), Rg(s.n, 1 + s.nc);
  NumericVector phi(s.n);
  for (int i = 0; i < s.n; ++i) {
    int yi = y[i];
    if (yi < 0 || yi >= J || !avail(i, yi)) { yi = 0; while (yi < J && !avail(i, yi)) ++yi; if (yi == J) yi = 0; }
    eval_case(i, s, theta, X, W, Z, off, avail, yi, e);
    for (int j = 0; j < J; ++j) { Pstar(i, j) = e.Pstar[j]; P(i, j) = e.P[j]; Vm(i, j) = e.V[j]; }
    Rg(i, 0) = e.pi0; for (int t = 0; t < s.nc; ++t) Rg(i, 1 + t) = e.pi[t];
    phi[i] = e.pi0;
  }
  return List::create(_["prob"] = Pstar, _["prob_outcome"] = P, _["phi"] = phi, _["regime"] = Rg, _["V"] = Vm);
}
