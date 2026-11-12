// mnl_ll.cpp -- the likelihood engine shared by mnl() and imnl(): the multinomial /
// conditional logit with case-specific and alternative-specific covariates,
// varying choice sets, offsets, and an optional inflated category.
//
// Model (case i = 0..n-1, alternatives j = 0..J-1, base b, availability a_ij):
//   V_ij  = off_ij + x_i' beta_j + w_ij' gamma,   beta_b = 0
//   P_ij  = a_ij exp(V_ij) / sum_k a_ik exp(V_ik)
//   split (k >= 0 only): eta_i = z_i' delta, phi_i = F(eta_i), F logistic or normal
//   P*_ij = phi_i P_ij + 1{j = k} (1 - phi_i)        (k < 0: P*_ij = P_ij)
//   l_i   = log P*_{i, y_i}
//
// Parameter vector theta, in this order:
//   beta   (p * (J-1))  alternative-major: for each non-base alternative in level
//                       order, the p case-specific coefficients
//   gamma  (q)          generic coefficients on the alternative-specific covariates
//   delta  (r)          split-equation coefficients (only when k >= 0)
//
// Designs: X is n x p (case-specific, intercept included when wanted); W is
// (n*J) x q with row i*J + j holding w_ij; Z is n x r; off and avail are n x J.
// y holds the chosen alternative (0-based) and must be available.
//
// Derivatives (exact). With d_j = 1{j = y} - P_j, P* = P*_{i,k}:
//   c  = 1                       (y != k)     c  = phi P_k / P*          (y == k)
//   dl/dV_j   = c d_j
//   u  = f/phi                   (y != k)     u  = -f (1 - P_k) / P*     (y == k)
//   dl/deta   = u
//   d2l/dV_j dV_m = dc/dV_m d_j - c P_j (1{j=m} - P_m),
//        dc/dV_m = 0 (y != k);  phi (1-phi) P_k (1{k=m} - P_m) / P*^2  (y == k)
//   d2l/dV_j deta = dc/deta d_j,  dc/deta = 0 (y != k);  f P_k / P*^2  (y == k)
//   d2l/deta^2    = f'/phi - u^2            (y != k)
//                 = -(1 - P_k) f'/P* - u^2  (y == k)
//   with f' = f (1 - 2 phi) (logit), -eta f (probit).
// The chain rule to theta uses dV_j/dbeta_m = x 1{j = m}, dV_j/dgamma = w_j,
// deta/ddelta = z. 1 - P_k is accumulated as the sum of the other available
// probabilities, log P*_k is a two-term log-sum-exp, and f/phi is formed on the
// log scale, so the extreme-split paths stay finite. Gradient, scores, and
// Hessian are validated against numDeriv in tests/testthat/test-derivatives.R.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
using namespace Rcpp;

namespace {

struct Spec {
  int n, J, p, q, r, base, k, link, K;
  std::vector<int> bpos;     // bpos[j]: offset of alternative j's beta block, -1 for the base
  int gpos, dpos;            // offsets of the gamma and delta blocks
};

Spec make_spec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
               const NumericMatrix& Z, const IntegerVector& y, int J, int base, int k, int link) {
  Spec s;
  s.n = y.size(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.r = (k >= 0) ? Z.ncol() : 0;
  s.base = base; s.k = k; s.link = link;
  s.bpos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  s.dpos = pos; pos += s.r;
  s.K = pos;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (X.nrow() != s.n || Z.nrow() != s.n) stop("internal: design rows do not match the response.");
  if (W.nrow() != s.n * J) stop("internal: W must have n * J rows.");
  return s;
}

inline double logF(double x, int link, bool lower) {
  return link == 1 ? R::plogis(x, 0.0, 1.0, lower ? 1 : 0, 1) : R::pnorm(x, 0.0, 1.0, lower ? 1 : 0, 1);
}
inline double fdens(double x, int link) {
  return link == 1 ? R::dlogis(x, 0.0, 1.0, 0) : R::dnorm(x, 0.0, 1.0, 0);
}
inline double lse2(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  double m = std::max(a, b);
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// Everything the entry points need for one case.
struct CaseEval {
  std::vector<double> P;      // outcome-stage probabilities (0 where unavailable)
  std::vector<double> Pstar;  // observed-category probabilities
  std::vector<double> dcdV;   // dc/dV_m (second order)
  double phi, logphi, log1mphi, f, eta;
  double lli;                 // log P*_{i, y_i}
  double c, u;                // first-order multipliers
  double dcdeta, dudeta;      // second-order multipliers
  void init(int J) { P.assign(J, 0.0); Pstar.assign(J, 0.0); dcdV.assign(J, 0.0); }
};

void eval_case(int i, const Spec& s, const NumericVector& theta, const NumericMatrix& X,
               const NumericMatrix& W, const NumericMatrix& Z, const NumericMatrix& off,
               const IntegerMatrix& avail, int yi, CaseEval& e, std::vector<double>& V) {
  const int J = s.J;
  double vmax = R_NegInf;
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) { V[j] = R_NegInf; continue; }
    double v = off(i, j);
    if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) v += X(i, c) * theta[s.bpos[j] + c];
    for (int c = 0; c < s.q; ++c) v += W(i * J + j, c) * theta[s.gpos + c];
    V[j] = v;
    if (v > vmax) vmax = v;
  }
  double denom = 0.0;
  for (int j = 0; j < J; ++j) if (V[j] != R_NegInf) denom += std::exp(V[j] - vmax);
  const double logdenom = std::log(denom);
  for (int j = 0; j < J; ++j) e.P[j] = (V[j] == R_NegInf) ? 0.0 : std::exp(V[j] - vmax - logdenom);
  const double logPy = V[yi] - vmax - logdenom;
  std::fill(e.dcdV.begin(), e.dcdV.end(), 0.0);
  e.dcdeta = 0.0; e.dudeta = 0.0;
  if (s.k < 0) {
    e.phi = 1.0; e.logphi = 0.0; e.log1mphi = R_NegInf; e.f = 0.0; e.eta = 0.0;
    e.lli = logPy; e.c = 1.0; e.u = 0.0;
    for (int j = 0; j < J; ++j) e.Pstar[j] = e.P[j];
    return;
  }
  double eta = 0.0;
  for (int c = 0; c < s.r; ++c) eta += Z(i, c) * theta[s.dpos + c];
  e.eta = eta;
  e.logphi = logF(eta, s.link, true);
  e.log1mphi = logF(eta, s.link, false);
  e.phi = std::exp(e.logphi);
  e.f = fdens(eta, s.link);
  const double logf = std::log(e.f);
  const double fprime = (s.link == 1) ? e.f * (1.0 - 2.0 * e.phi) : -eta * e.f;
  if (!avail(i, s.k)) {
    // the inflated alternative is not in this case's choice set: regime 0
    // ("always k") is impossible, the case is in the choice regime with
    // certainty, its contribution is the MNL one, and it carries no
    // information about delta (all split-equation derivatives vanish)
    e.phi = 1.0; e.logphi = 0.0; e.log1mphi = R_NegInf; e.f = 0.0;
    e.lli = logPy; e.c = 1.0; e.u = 0.0;
    for (int j = 0; j < J; ++j) e.Pstar[j] = e.P[j];
    return;
  }
  const double logPk = V[s.k] - vmax - logdenom;
  const double logPstar_k = lse2(e.log1mphi, e.logphi + logPk);
  for (int j = 0; j < J; ++j) e.Pstar[j] = (j == s.k) ? std::exp(logPstar_k) : e.phi * e.P[j];
  if (yi != s.k) {
    e.lli = e.logphi + logPy;
    e.c = 1.0;
    e.u = std::exp(logf - e.logphi);                       // f / phi
    const double fp_over_phi = (s.link == 1) ? e.u * (1.0 - 2.0 * e.phi) : -eta * e.u;
    e.dudeta = fp_over_phi - e.u * e.u;
  } else {
    e.lli = logPstar_k;
    const double Pstar = std::exp(logPstar_k);
    e.c = std::exp(e.logphi + logPk - logPstar_k);          // phi P_k / P*, in [0, 1]
    double one_minus_Pk = 0.0;
    for (int j = 0; j < J; ++j) if (j != s.k) one_minus_Pk += e.P[j];
    e.u = -e.f * one_minus_Pk / Pstar;
    const double Pk = e.P[s.k];
    const double phi1mphi = std::exp(e.logphi + e.log1mphi);
    for (int m = 0; m < J; ++m) if (V[m] != R_NegInf)
      e.dcdV[m] = phi1mphi * Pk * ((m == s.k ? 1.0 : 0.0) - e.P[m]) / (Pstar * Pstar);
    e.dcdeta = e.f * Pk / (Pstar * Pstar);
    e.dudeta = -one_minus_Pk * fprime / Pstar - e.u * e.u;
  }
}

// score of case i into out[0..K-1] (overwritten)
void score_case(int i, const Spec& s, const CaseEval& e, const NumericMatrix& X, const NumericMatrix& W,
                const NumericMatrix& Z, const IntegerMatrix& avail, int yi, std::vector<double>& out) {
  std::fill(out.begin(), out.end(), 0.0);
  const int J = s.J;
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    const double r = e.c * ((j == yi ? 1.0 : 0.0) - e.P[j]);
    if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) out[s.bpos[j] + c] += r * X(i, c);
    for (int c = 0; c < s.q; ++c) out[s.gpos + c] += r * W(i * J + j, c);
  }
  if (s.k >= 0) for (int c = 0; c < s.r; ++c) out[s.dpos + c] = e.u * Z(i, c);
}

} // namespace

// Per-case log-likelihood contributions log P*_{i, y_i} (unweighted).
// [[Rcpp::export]]
NumericVector pc_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                         NumericMatrix off, IntegerMatrix avail, IntegerVector y,
                         int J, int base, int k, int link) {
  Spec s = make_spec(theta, X, W, Z, y, J, base, k, link);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  CaseEval e; e.init(J);
  std::vector<double> V(J);
  NumericVector out(s.n);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, Z, off, avail, y[i], e, V);
    out[i] = e.lli;
  }
  return out;
}

// Gradient of the weighted log-likelihood sum_i w_i log P*_{i, y_i}.
// [[Rcpp::export]]
NumericVector pc_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                          NumericMatrix off, IntegerMatrix avail, IntegerVector y, NumericVector w,
                          int J, int base, int k, int link) {
  Spec s = make_spec(theta, X, W, Z, y, J, base, k, link);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  if (w.size() != s.n) stop("internal: weights do not match the cases.");
  CaseEval e; e.init(J);
  std::vector<double> V(J), sc(s.K);
  NumericVector g(s.K);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, Z, off, avail, y[i], e, V);
    score_case(i, s, e, X, W, Z, avail, y[i], sc);
    for (int c = 0; c < s.K; ++c) g[c] += w[i] * sc[c];
  }
  return g;
}

// Per-case scores s_i (unweighted), n x K.
// [[Rcpp::export]]
NumericMatrix pc_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                            NumericMatrix off, IntegerMatrix avail, IntegerVector y,
                            int J, int base, int k, int link) {
  Spec s = make_spec(theta, X, W, Z, y, J, base, k, link);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  CaseEval e; e.init(J);
  std::vector<double> V(J), sc(s.K);
  NumericMatrix S(s.n, s.K);
  for (int i = 0; i < s.n; ++i) {
    eval_case(i, s, theta, X, W, Z, off, avail, y[i], e, V);
    score_case(i, s, e, X, W, Z, avail, y[i], sc);
    for (int c = 0; c < s.K; ++c) S(i, c) = sc[c];
  }
  return S;
}

// Hessian of the weighted log-likelihood (K x K, negative semi-definite at the
// maximum). The per-alternative design rows D_j = dV_j/dtheta are handled as
// sparse vectors, so unit fixed-effect dummies cost only their nonzeros.
// [[Rcpp::export]]
NumericMatrix pc_hess_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                          NumericMatrix off, IntegerMatrix avail, IntegerVector y, NumericVector w,
                          int J, int base, int k, int link) {
  Spec s = make_spec(theta, X, W, Z, y, J, base, k, link);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  CaseEval e; e.init(J);
  std::vector<double> V(J), d(J);
  std::vector< std::vector< std::pair<int, double> > > D(J);
  NumericMatrix H(s.K, s.K);
  for (int i = 0; i < s.n; ++i) {
    const int yi = y[i];
    eval_case(i, s, theta, X, W, Z, off, avail, yi, e, V);
    const double wi = w[i];
    for (int j = 0; j < J; ++j) {
      D[j].clear();
      if (!avail(i, j)) continue;
      d[j] = (j == yi ? 1.0 : 0.0) - e.P[j];
      if (s.bpos[j] >= 0) for (int c = 0; c < s.p; ++c) { const double v = X(i, c); if (v != 0.0) D[j].push_back(std::make_pair(s.bpos[j] + c, v)); }
      for (int c = 0; c < s.q; ++c) { const double v = W(i * J + j, c); if (v != 0.0) D[j].push_back(std::make_pair(s.gpos + c, v)); }
    }
    // V-V block: d2l/dV_j dV_m = dc/dV_m d_j - c P_j (1{j=m} - P_m)
    for (int j = 0; j < J; ++j) {
      if (!avail(i, j) || D[j].empty()) continue;
      for (int m = 0; m < J; ++m) {
        if (!avail(i, m) || D[m].empty()) continue;
        const double h = wi * (e.dcdV[m] * d[j] - e.c * e.P[j] * ((j == m ? 1.0 : 0.0) - e.P[m]));
        if (h == 0.0) continue;
        for (size_t a = 0; a < D[j].size(); ++a)
          for (size_t b = 0; b < D[m].size(); ++b)
            H(D[j][a].first, D[m][b].first) += h * D[j][a].second * D[m][b].second;
      }
    }
    if (s.k >= 0 && s.r > 0) {
      // V-eta cross terms and the eta-eta block
      for (int j = 0; j < J; ++j) {
        if (!avail(i, j) || D[j].empty()) continue;
        const double hv = wi * e.dcdeta * d[j];
        if (hv == 0.0) continue;
        for (size_t a = 0; a < D[j].size(); ++a)
          for (int c = 0; c < s.r; ++c) {
            const double v = hv * D[j][a].second * Z(i, c);
            H(D[j][a].first, s.dpos + c) += v;
            H(s.dpos + c, D[j][a].first) += v;
          }
      }
      for (int c = 0; c < s.r; ++c)
        for (int c2 = 0; c2 < s.r; ++c2)
          H(s.dpos + c, s.dpos + c2) += wi * e.dudeta * Z(i, c) * Z(i, c2);
    }
  }
  return H;
}

// Probabilities: P* (observed categories), P (outcome stage), phi (P(s_i = 1)),
// and the linear indices V (n x J, -Inf where unavailable). `y` is only used
// for its length here; pass any valid alternative.
// [[Rcpp::export]]
List pc_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix Z,
                  NumericMatrix off, IntegerMatrix avail, IntegerVector y,
                  int J, int base, int k, int link) {
  Spec s = make_spec(theta, X, W, Z, y, J, base, k, link);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  CaseEval e; e.init(J);
  std::vector<double> V(J);
  NumericMatrix Pstar(s.n, J), P(s.n, J), Vm(s.n, J);
  NumericVector phi(s.n);
  for (int i = 0; i < s.n; ++i) {
    int yi = y[i];
    if (yi < 0 || yi >= J || !avail(i, yi)) { yi = 0; while (yi < J && !avail(i, yi)) ++yi; if (yi == J) yi = 0; }
    eval_case(i, s, theta, X, W, Z, off, avail, yi, e, V);
    for (int j = 0; j < J; ++j) { Pstar(i, j) = e.Pstar[j]; P(i, j) = e.P[j]; Vm(i, j) = V[j]; }
    phi[i] = e.phi;
  }
  return List::create(_["prob"] = Pstar, _["prob_outcome"] = P, _["phi"] = phi, _["V"] = Vm);
}
