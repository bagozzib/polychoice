// remnl_ll.cpp -- the multinomial logit with random alternative-specific
// intercepts (a mixed logit whose only random coefficients are the constants):
//
//   V_ij(u) = V_ij + u_{g(i) j},   V_ij as in mnl_ll.cpp,
//   u_g = (u_{g j})_{j != base} ~ N(0, diag(sigma_j^2)), independent across
//   units g and alternatives, u_{g base} = 0.
//   L_g = int prod_{i in g} P_i(y_i | u) phi_sigma(u) du,   P_i the MNL
//         probability at the shifted utilities.
//
// Adaptive Gauss-Hermite quadrature (Liu and Pierce 1994; Pinheiro and Bates
// 1995). The unit log-integrand h_g(u) = sum_i log P_i(y_i | u) + log phi_sigma(u)
// is strictly concave (the MNL log-likelihood is concave in the utilities and
// the prior is log-concave), so its mode u_g^ is found by damped Newton steps
// from u = 0, with the m x m Hessian H_g. With -H_g = L L' and C = L^{-T}
// (so C C' = (-H_g)^{-1}), the nodes are u_q = u_g^ + C z_q for the
// probabilists' product Gauss-Hermite points z_q (Q x m, weights w_q), and
//   L_g = int exp(h_g(u)) du = |C| sum_q w_q exp(h_g(u_q) + z_q'z_q / 2) (2 pi)^{m/2},
// a rule that is exact when the integrand is Gaussian and whose error shrinks
// with the number of cases per unit instead of growing (a fixed rule centred at
// zero with nodes sigma z_q aliases once the integrand is narrower than the node
// spacing: the ledger records -22 log-likelihood units and -25 percent in sigma
// at 50 cases per unit with 12 fixed nodes against lme4::glmer).
//
// The observation is the unit: pc_re_lli_cpp returns one log-likelihood per
// unit, pc_re_scores_cpp one score row per unit, the mixture-weighted average of
// the per-node derivatives with the nodes held fixed (the integral does not
// depend on them; the dependence of u_q on theta through the mode is of the
// order of the quadrature error). Parameter vector theta: beta (p (J-1)), gamma
// (q), then log sigma_j for the non-base alternatives in level order;
//   d h / d V_ij = 1{y_i = j} - P_ij(u_q)  (beta and gamma scores as in the MNL),
//   d h / d log sigma_j = u_qj^2 / sigma_j^2 - 1   (the prior).
// pc_re_probs_cpp returns the population-averaged probabilities
// sum_q w_q P_ij(sigma z_q) of every case on the rule it is given (no data
// narrow this integrand, so the R side passes a finer rule than the
// likelihood's: .pc_re_pred_rule) and the utilities V_ij. Cases are passed grouped
// by unit: uidx holds the case indices ordered by unit, ustart[g]..ustart[g+1]-1
// the positions of unit g in uidx.
// Validated against numDeriv and, at J = 2, against lme4::glmer in
// tests/testthat/test-re.R.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

namespace {

struct RESpec {
  int n, J, p, q, m, base, K, G, Q;
  std::vector<int> bpos, spos, npos;   // beta block start per alternative, log-sigma index, non-base position (column of Z)
  int gpos;
};

RESpec re_spec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W, const NumericMatrix& off,
               const IntegerMatrix& avail, const IntegerVector& y, int J, int base,
               const IntegerVector& ustart, const NumericMatrix& Z, const NumericVector& logw) {
  RESpec s;
  s.n = y.size(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.m = J - 1; s.base = base;
  s.G = ustart.size() - 1; s.Q = Z.nrow();
  s.bpos.assign(J, -1); s.spos.assign(J, -1); s.npos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  int np = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.spos[j] = pos + np; s.npos[j] = np; ++np; }
  pos += s.m;
  s.K = pos;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (X.nrow() != s.n || W.nrow() != s.n * J) stop("internal: design rows do not match the response.");
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  if (Z.ncol() != s.m || logw.size() != s.Q) stop("internal: quadrature nodes do not match the alternatives.");
  return s;
}

// base utilities (without u) of one case into V (length J; -Inf where unavailable)
inline void re_utilities(int i, const RESpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
                         const NumericMatrix& off, const IntegerMatrix& avail, std::vector<double>& V) {
  for (int j = 0; j < s.J; ++j) {
    if (!avail(i, j)) { V[j] = R_NegInf; continue; }
    double v = off(i, j);
    if (s.bpos[j] >= 0) for (int cc = 0; cc < s.p; ++cc) v += X(i, cc) * theta[s.bpos[j] + cc];
    for (int cc = 0; cc < s.q; ++cc) v += W(i * s.J + j, cc) * theta[s.gpos + cc];
    V[j] = v;
  }
}

// MNL probabilities of one case at the shifted utilities; returns log P(y) and fills P
inline double re_case_probs(const RESpec& s, const std::vector<double>& Vb, const IntegerMatrix& avail, int i, int yi,
                            const std::vector<double>& u, std::vector<double>& P) {
  const int J = s.J;
  double mx = R_NegInf;
  std::vector<double>& Vu = P;                     // reuse storage
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) { Vu[j] = R_NegInf; continue; }
    Vu[j] = Vb[j] + (s.npos[j] >= 0 ? u[s.npos[j]] : 0.0);
    if (Vu[j] > mx) mx = Vu[j];
  }
  double se = 0.0;
  for (int j = 0; j < J; ++j) if (avail(i, j)) se += std::exp(Vu[j] - mx);
  const double lse = mx + std::log(se);
  const double lpy = Vu[yi] - lse;
  for (int j = 0; j < J; ++j) P[j] = avail(i, j) ? std::exp(Vu[j] - lse) : 0.0;
  return lpy;
}

// Cholesky of a small symmetric positive definite matrix A (m x m in a flat
// vector, symmetric, so either order); Lo lower triangular, row-major (index
// i*m + k); false if it fails
inline bool chol_small(int m, const std::vector<double>& A, std::vector<double>& Lo) {
  Lo.assign((size_t) m * m, 0.0);
  for (int j = 0; j < m; ++j) {
    double d = A[j * m + j];
    for (int k = 0; k < j; ++k) d -= Lo[j * m + k] * Lo[j * m + k];
    if (!(d > 0.0)) return false;
    Lo[j * m + j] = std::sqrt(d);
    for (int i = j + 1; i < m; ++i) {
      double v = A[i * m + j];
      for (int k = 0; k < j; ++k) v -= Lo[i * m + k] * Lo[j * m + k];
      Lo[i * m + j] = v / Lo[j * m + j];
    }
  }
  return true;
}
// solve L L' x = b in place (L lower triangular, row-major index i*m + k)
inline void chol_solve(int m, const std::vector<double>& Lo, std::vector<double>& b) {
  for (int i = 0; i < m; ++i) { double v = b[i]; for (int k = 0; k < i; ++k) v -= Lo[i * m + k] * b[k]; b[i] = v / Lo[i * m + i]; }
  for (int i = m - 1; i >= 0; --i) { double v = b[i]; for (int k = i + 1; k < m; ++k) v -= Lo[k * m + i] * b[k]; b[i] = v / Lo[i * m + i]; }
}

// log-likelihood of unit g and, when asked, its score
double re_unit(int g, const RESpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
               const NumericMatrix& off, const IntegerMatrix& avail, const IntegerVector& y, const IntegerVector& uidx,
               const IntegerVector& ustart, const NumericMatrix& Z, const NumericVector& logw, const std::vector<double>& sigma,
               bool want_grad, std::vector<double>& score) {
  const int J = s.J, Q = s.Q, m = s.m;
  const int from = ustart[g], to = ustart[g + 1];
  const int T = to - from;
  std::vector< std::vector<double> > Vb(T, std::vector<double>(J));
  for (int t = 0; t < T; ++t) re_utilities(uidx[from + t], s, theta, X, W, off, avail, Vb[t]);
  std::vector<double> P(J), u(m, 0.0), gvec(m), Hm((size_t) m * m), Lo, step(m), unew(m);
  // the unit log-integrand h(u) (without constants) and its gradient / Hessian
  auto h_of = [&](const std::vector<double>& uu) {
    double h = 0.0;
    for (int t = 0; t < T; ++t) h += re_case_probs(s, Vb[t], avail, uidx[from + t], y[uidx[from + t]], uu, P);
    for (int k = 0; k < m; ++k) h -= 0.5 * uu[k] * uu[k] / (sigma[k] * sigma[k]);
    return h;
  };
  // Newton for the mode from u = 0 (concave: damped steps; about one unit of u
  // per step in the logistic tail, so the cap of 100 is never reached in practice)
  double h = h_of(u);
  for (int it = 0; it < 100; ++it) {
    std::fill(gvec.begin(), gvec.end(), 0.0); std::fill(Hm.begin(), Hm.end(), 0.0);
    for (int t = 0; t < T; ++t) {
      const int i = uidx[from + t], yi = y[i];
      re_case_probs(s, Vb[t], avail, i, yi, u, P);
      for (int j = 0; j < J; ++j) {
        if (s.npos[j] < 0 || !avail(i, j)) continue;
        const int a = s.npos[j];
        gvec[a] += (j == yi ? 1.0 : 0.0) - P[j];
        for (int k = 0; k < J; ++k) {
          if (s.npos[k] < 0 || !avail(i, k)) continue;
          const int b = s.npos[k];
          Hm[a * m + b] += P[j] * ((j == k ? 1.0 : 0.0) - P[k]);        // minus the Hessian of the log-likelihood part
        }
      }
    }
    for (int k = 0; k < m; ++k) { gvec[k] -= u[k] / (sigma[k] * sigma[k]); Hm[k * m + k] += 1.0 / (sigma[k] * sigma[k]); }
    double gmax = 0.0; for (int k = 0; k < m; ++k) gmax = std::max(gmax, std::fabs(gvec[k]));
    if (gmax < 1e-9) break;
    if (!chol_small(m, Hm, Lo)) break;                 // cannot happen for a strictly concave h; be safe
    step = gvec; chol_solve(m, Lo, step);              // Newton direction (-H)^{-1} g
    double lam = 1.0; bool ok = false;
    for (int hv = 0; hv < 20; ++hv) {
      for (int k = 0; k < m; ++k) unew[k] = u[k] + lam * step[k];
      const double hn = h_of(unew);
      if (std::isfinite(hn) && hn >= h - 1e-12 * std::fabs(h)) { u = unew; h = hn; ok = true; break; }
      lam *= 0.5;
    }
    if (!ok) break;
  }
  // curvature at the mode: -H = L L', C = L^{-T}, |C| = 1 / prod L_kk
  std::fill(Hm.begin(), Hm.end(), 0.0);
  for (int t = 0; t < T; ++t) {
    const int i = uidx[from + t];
    re_case_probs(s, Vb[t], avail, i, y[i], u, P);
    for (int j = 0; j < J; ++j) {
      if (s.npos[j] < 0 || !avail(i, j)) continue;
      const int a = s.npos[j];
      for (int k = 0; k < J; ++k) {
        if (s.npos[k] < 0 || !avail(i, k)) continue;
        Hm[a * m + s.npos[k]] += P[j] * ((j == k ? 1.0 : 0.0) - P[k]);
      }
    }
  }
  for (int k = 0; k < m; ++k) Hm[k * m + k] += 1.0 / (sigma[k] * sigma[k]);
  if (!chol_small(m, Hm, Lo)) { Lo.assign((size_t) m * m, 0.0); for (int k = 0; k < m; ++k) Lo[k * m + k] = 1.0 / sigma[k]; }
  double logdetC = 0.0;
  for (int k = 0; k < m; ++k) logdetC -= std::log(Lo[k * m + k]);
  // the quadrature: u_q = u^ + C z_q with C = L^{-T}, i.e. solve L' c = z
  std::vector<double> lp(Q), uq(m), c(m);
  std::vector<double> sc;
  if (want_grad) sc.assign((size_t) Q * s.K, 0.0);
  const double half_log_2pi_m = 0.5 * m * std::log(2.0 * M_PI);
  for (int q = 0; q < Q; ++q) {
    for (int k = 0; k < m; ++k) c[k] = Z(q, k);
    for (int i = m - 1; i >= 0; --i) { double v = c[i]; for (int k = i + 1; k < m; ++k) v -= Lo[k * m + i] * c[k]; c[i] = v / Lo[i * m + i]; }
    double zz = 0.0, lprior = 0.0;
    for (int k = 0; k < m; ++k) {
      uq[k] = u[k] + c[k];
      zz += Z(q, k) * Z(q, k);
      lprior += -std::log(sigma[k]) - 0.5 * std::log(2.0 * M_PI) - 0.5 * uq[k] * uq[k] / (sigma[k] * sigma[k]);
    }
    double acc = logw[q] + logdetC + half_log_2pi_m + 0.5 * zz + lprior;
    double* sq = want_grad ? &sc[(size_t) q * s.K] : 0;
    for (int t = 0; t < T; ++t) {
      const int i = uidx[from + t], yi = y[i];
      acc += re_case_probs(s, Vb[t], avail, i, yi, uq, P);
      if (want_grad) {
        for (int k = 0; k < J; ++k) {
          if (!avail(i, k)) continue;
          const double d = (k == yi ? 1.0 : 0.0) - P[k];
          if (s.bpos[k] >= 0) for (int cc = 0; cc < s.p; ++cc) sq[s.bpos[k] + cc] += d * X(i, cc);
          for (int cc = 0; cc < s.q; ++cc) sq[s.gpos + cc] += d * W(i * J + k, cc);
        }
      }
    }
    if (want_grad) for (int j = 0; j < J; ++j) if (s.npos[j] >= 0) { const int a = s.npos[j]; sq[s.spos[j]] += uq[a] * uq[a] / (sigma[a] * sigma[a]) - 1.0; }
    lp[q] = acc;
  }
  double mx = R_NegInf;
  for (int q = 0; q < Q; ++q) if (lp[q] > mx) mx = lp[q];
  double se = 0.0;
  for (int q = 0; q < Q; ++q) se += std::exp(lp[q] - mx);
  const double ll = mx + std::log(se);
  if (want_grad) {
    std::fill(score.begin(), score.end(), 0.0);
    for (int q = 0; q < Q; ++q) {
      const double pi = std::exp(lp[q] - ll);
      if (pi == 0.0) continue;
      const double* sq = &sc[(size_t) q * s.K];
      for (int k = 0; k < s.K; ++k) score[k] += pi * sq[k];
    }
  }
  return ll;
}

std::vector<double> re_sigma(const RESpec& s, const NumericVector& theta) {
  std::vector<double> sigma(s.m);
  for (int j = 0; j < s.J; ++j) if (s.spos[j] >= 0) sigma[s.npos[j]] = std::exp(theta[s.spos[j]]);
  return sigma;
}

} // namespace

// [[Rcpp::export]]
NumericVector pc_re_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                            IntegerVector y, int J, int base, IntegerVector uidx, IntegerVector ustart, NumericMatrix Z, NumericVector logw) {
  RESpec s = re_spec(theta, X, W, off, avail, y, J, base, ustart, Z, logw);
  std::vector<double> sigma = re_sigma(s, theta), score;
  NumericVector out(s.G);
  for (int g = 0; g < s.G; ++g) {
    if (g % 500 == 0) Rcpp::checkUserInterrupt();
    out[g] = re_unit(g, s, theta, X, W, off, avail, y, uidx, ustart, Z, logw, sigma, false, score);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pc_re_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                             IntegerVector y, NumericVector w, int J, int base, IntegerVector uidx, IntegerVector ustart, NumericMatrix Z, NumericVector logw) {
  RESpec s = re_spec(theta, X, W, off, avail, y, J, base, ustart, Z, logw);
  if (w.size() != s.G) stop("internal: unit weights do not match the units.");
  std::vector<double> sigma = re_sigma(s, theta), score(s.K);
  NumericVector g(s.K);
  for (int u = 0; u < s.G; ++u) {
    if (u % 500 == 0) Rcpp::checkUserInterrupt();
    re_unit(u, s, theta, X, W, off, avail, y, uidx, ustart, Z, logw, sigma, true, score);
    for (int k = 0; k < s.K; ++k) g[k] += w[u] * score[k];
  }
  return g;
}

// [[Rcpp::export]]
NumericMatrix pc_re_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                               IntegerVector y, int J, int base, IntegerVector uidx, IntegerVector ustart, NumericMatrix Z, NumericVector logw) {
  RESpec s = re_spec(theta, X, W, off, avail, y, J, base, ustart, Z, logw);
  std::vector<double> sigma = re_sigma(s, theta), score(s.K);
  NumericMatrix S(s.G, s.K);
  for (int u = 0; u < s.G; ++u) {
    if (u % 500 == 0) Rcpp::checkUserInterrupt();
    re_unit(u, s, theta, X, W, off, avail, y, uidx, ustart, Z, logw, sigma, true, score);
    for (int k = 0; k < s.K; ++k) S(u, k) = score[k];
  }
  return S;
}

// population-averaged probabilities of every case (prior rule) and the utilities without u
// [[Rcpp::export]]
List pc_re_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                     IntegerVector y, int J, int base, IntegerVector uidx, IntegerVector ustart, NumericMatrix Z, NumericVector logw) {
  RESpec s = re_spec(theta, X, W, off, avail, y, J, base, ustart, Z, logw);
  std::vector<double> sigma = re_sigma(s, theta);
  NumericMatrix P(s.n, J), Vm(s.n, J);
  std::vector<double> V(J), Vu(J);
  for (int i = 0; i < s.n; ++i) {
    if (i % 2000 == 0) Rcpp::checkUserInterrupt();
    re_utilities(i, s, theta, X, W, off, avail, V);
    for (int j = 0; j < J; ++j) { Vm(i, j) = V[j]; P(i, j) = 0.0; }
    for (int q = 0; q < s.Q; ++q) {
      double mx = R_NegInf;
      for (int j = 0; j < J; ++j) {
        if (!avail(i, j)) { Vu[j] = R_NegInf; continue; }
        Vu[j] = V[j] + (s.npos[j] >= 0 ? sigma[s.npos[j]] * Z(q, s.npos[j]) : 0.0);
        if (Vu[j] > mx) mx = Vu[j];
      }
      double se = 0.0;
      for (int j = 0; j < J; ++j) if (avail(i, j)) se += std::exp(Vu[j] - mx);
      const double wq = std::exp(logw[q]);
      for (int j = 0; j < J; ++j) if (avail(i, j)) P(i, j) += wq * std::exp(Vu[j] - mx) / se;
    }
  }
  return List::create(_["prob"] = P, _["V"] = Vm);
}
