// mnp_ll.cpp -- the multinomial probit: utilities U_ij = V_ij + eps_ij with
// eps_i multivariate normal. Only utility differences are identified, so the
// model is parameterized through the covariance of the differences from the
// base alternative, eta_k = eps_k - eps_b (k != b), Sigma_tilde = L L' with L
// lower triangular, L[0][0] = 1 (scale normalization), the other diagonal
// elements positive (internally on the log scale), off-diagonals free:
// (J-1)J/2 - 1 free parameters, as in Train (2009, section 5.2) and mlogit.
//
//   P(y = j) = P(eps_k - eps_j <= V_j - V_k for all available k != j)
//            = Phi_m(d_j; Omega_j),  d_jk = V_j - V_k,  Omega_j = B B',  B = A_j L,
// where A_j maps the differences from the base to the differences from j and
// m is the number of available competitors. m = 1: Phi; m = 2 (three
// alternatives): the exact bivariate normal of pc_bvn.h; m >= 3: the GHK
// simulator (Geweke-Hajivassiliou-Keane; Train 2009, section 5.6) with R
// Halton draws per case (bases 2, 3, 5, ..., the first 100 points skipped;
// each case uses its own block of points, held fixed across evaluations so
// the simulated likelihood is smooth).
//
// Parameter vector theta: beta (p (J-1)), gamma (q), then the free elements of
// L in row-major lower-triangular order (log for the diagonal).
//
// Derivatives are exact (given the draws), by forward-mode differentiation in
// mc + nc directions (the utility differences d_t and the structural
// parameters c_t): dOmega = dB B' + B dB', the Cholesky differential
// dC = C Phi(C^{-1} dOmega C^{-T}) (Phi: lower triangle with halved diagonal),
// and through the GHK recursion b_s = (d_s - sum_{l<s} C_sl e_l) / C_ss,
// e_s = Phi^{-1}(u_s Phi(b_s)):
//   db_s = (dd_s - sum_{l<s} (dC_sl e_l + C_sl de_l)) / C_ss - b_s dC_ss / C_ss,
//   de_s = u_s phi(b_s) / phi(e_s) db_s,   d log P_r = sum_s phi(b_s)/Phi(b_s) db_s.
// Validated against numDeriv in tests/testthat/test-mnp.R.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "pc_bvn.h"
using namespace Rcpp;

namespace {

struct PSpec {
  int n, J, p, q, base, m, nc, K, R;
  std::vector<int> bpos;
  int gpos, cpos;
  std::vector<int> pos;              // pos[k]: position of alternative k among the non-base alternatives, -1 for the base
  std::vector< std::pair<int,int> > cidx;   // free L elements (row, col) in theta order
};

PSpec mnp_spec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W, const IntegerVector& y, int J, int base, int R) {
  PSpec s;
  s.n = y.size(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.base = base; s.m = J - 1; s.R = R;
  s.bpos.assign(J, -1); s.pos.assign(J, -1);
  int pos = 0, np = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; s.pos[j] = np++; }
  s.gpos = pos; pos += s.q;
  s.cpos = pos;
  for (int r = 0; r < s.m; ++r) for (int c = 0; c <= r; ++c) if (!(r == 0 && c == 0)) s.cidx.push_back(std::make_pair(r, c));
  s.nc = s.cidx.size(); pos += s.nc;
  s.K = pos;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (X.nrow() != s.n || W.nrow() != s.n * J) stop("internal: design rows do not match the response.");
  if (J > 21) stop("internal: the GHK simulator supports at most 21 alternatives (20 Halton bases).");
  return s;
}

// L (row-major m x m) and the derivative of L with respect to each structural parameter
void build_L(const PSpec& s, const double* cpar, std::vector<double>& L, std::vector< std::vector<double> >& dL) {
  const int m = s.m;
  L.assign(m * m, 0.0); L[0] = 1.0;
  dL.assign(s.nc, std::vector<double>(m * m, 0.0));
  for (int t = 0; t < s.nc; ++t) {
    const int r = s.cidx[t].first, c = s.cidx[t].second;
    const double v = (r == c) ? std::exp(cpar[t]) : cpar[t];
    L[r * m + c] = v;
    dL[t][r * m + c] = (r == c) ? v : 1.0;       // d exp(c)/dc = exp(c)
  }
}

inline double radical_inverse(long long idx, int b) {
  double f = 1.0, r = 0.0;
  while (idx > 0) { f /= b; r += f * (idx % b); idx /= b; }
  return r;
}
static const int PRIMES[20] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

bool chol_small(int m, const std::vector<double>& A, std::vector<double>& C) {
  C.assign(m * m, 0.0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = A[i * m + j];
      for (int k = 0; k < j; ++k) s -= C[i * m + k] * C[j * m + k];
      if (i == j) { if (s <= 1e-300) return false; C[i * m + i] = std::sqrt(s); }
      else C[i * m + j] = s / C[j * m + j];
    }
  }
  return true;
}

// dC = C Phi(C^{-1} dOm C^{-T}) for a lower-triangular Cholesky factor C (row-major m x m)
void chol_diff(int m, const std::vector<double>& C, const std::vector<double>& dOm, std::vector<double>& dC) {
  // Y = C^{-1} dOm (forward substitution on each column), then M = Y C^{-T} = (C^{-1} Y')' ... compute M = C^{-1} dOm C^{-T}
  std::vector<double> Y(m * m, 0.0), M(m * m, 0.0);
  for (int col = 0; col < m; ++col)
    for (int i = 0; i < m; ++i) {
      double s = dOm[i * m + col];
      for (int k = 0; k < i; ++k) s -= C[i * m + k] * Y[k * m + col];
      Y[i * m + col] = s / C[i * m + i];
    }
  // M = Y C^{-T}: solve M C' = Y, i.e. for each row of Y: C z = y' -> z = M row
  for (int row = 0; row < m; ++row)
    for (int i = 0; i < m; ++i) {
      double s = Y[row * m + i];
      for (int k = 0; k < i; ++k) s -= C[i * m + k] * M[row * m + k];
      M[row * m + i] = s / C[i * m + i];
    }
  // Phi(M): lower triangle, halved diagonal
  std::vector<double> Ph(m * m, 0.0);
  for (int i = 0; i < m; ++i) for (int j = 0; j <= i; ++j) Ph[i * m + j] = (i == j) ? 0.5 * M[i * m + j] : M[i * m + j];
  dC.assign(m * m, 0.0);
  for (int i = 0; i < m; ++i) for (int j = 0; j <= i; ++j) {
    double s = 0.0; for (int k = 0; k < m; ++k) s += C[i * m + k] * Ph[k * m + j];
    dC[i * m + j] = s;
  }
}

struct PCase {
  std::vector<double> V, P, dV, dc;
  double lli;
};

// log P(j) and, when want_grad, d log P / d d_t (competitors t) and d log P / d c_t (structural)
double mnp_logprob(int i, int j, const PSpec& s, const std::vector<double>& V, const std::vector<double>& L,
                   const std::vector< std::vector<double> >& dL, const IntegerMatrix& avail, bool want_grad,
                   std::vector<int>& comp, std::vector<double>& dlogP_dd, std::vector<double>& dlogP_dc) {
  const int J = s.J, m = s.m, nc = s.nc;
  comp.clear();
  for (int k = 0; k < J; ++k) if (k != j && avail(i, k)) comp.push_back(k);
  const int mc = comp.size();
  dlogP_dd.assign(mc, 0.0); dlogP_dc.assign(nc, 0.0);
  if (mc == 0) return 0.0;
  // B = A_j L restricted to the competitors; dB_t = A_j dL_t
  auto build_B = [&](const std::vector<double>& LL, std::vector<double>& BB) {
    BB.assign(mc * m, 0.0);
    for (int t = 0; t < mc; ++t) {
      const int k = comp[t];
      for (int c = 0; c < m; ++c) {
        double v = 0.0;
        if (k != s.base) v += LL[s.pos[k] * m + c];
        if (j != s.base) v -= LL[s.pos[j] * m + c];
        BB[t * m + c] = v;
      }
    }
  };
  std::vector<double> B, d(mc); build_B(L, B);
  for (int t = 0; t < mc; ++t) d[t] = V[j] - V[comp[t]];
  std::vector<double> Om(mc * mc, 0.0);
  for (int a = 0; a < mc; ++a) for (int b = 0; b < mc; ++b) { double v = 0.0; for (int c = 0; c < m; ++c) v += B[a * m + c] * B[b * m + c]; Om[a * mc + b] = v; }
  std::vector< std::vector<double> > dOm;
  if (want_grad) {
    dOm.assign(nc, std::vector<double>(mc * mc, 0.0));
    std::vector<double> dB;
    for (int t = 0; t < nc; ++t) {
      build_B(dL[t], dB);
      for (int a = 0; a < mc; ++a) for (int b = 0; b < mc; ++b) {
        double v = 0.0; for (int c = 0; c < m; ++c) v += dB[a * m + c] * B[b * m + c] + B[a * m + c] * dB[b * m + c];
        dOm[t][a * mc + b] = v;
      }
    }
  }
  if (mc == 1) {
    const double sd = std::sqrt(Om[0]);
    const double z = d[0] / sd;
    const double lp = R::pnorm(z, 0.0, 1.0, 1, 1);
    if (want_grad) {
      const double lam = std::exp(R::dnorm(z, 0.0, 1.0, 1) - lp);
      dlogP_dd[0] = lam / sd;
      for (int t = 0; t < nc; ++t) { const double dsd = dOm[t][0] / (2.0 * sd); dlogP_dc[t] = lam * (-z * dsd / sd); }
    }
    return lp;
  }
  if (mc == 2) {
    const double s1 = std::sqrt(Om[0]), s2 = std::sqrt(Om[3]);
    double rho = Om[1] / (s1 * s2);
    rho = std::max(-0.999999, std::min(0.999999, rho));
    const double a = d[0] / s1, b = d[1] / s2;
    // a probability below the floor (an optimizer excursion far into a tail)
    // is held at the floor, where the log-likelihood is flat: its gradient is
    // 0 rather than a derivative divided by the floor
    const double Praw = pcbvn::pbvn(a, b, rho);
    const double P = std::max(Praw, 1e-300);
    const double lp = std::log(P);
    if (want_grad && !(Praw > 1e-300)) {
      dlogP_dd[0] = dlogP_dd[1] = 0.0;
      for (int t = 0; t < nc; ++t) dlogP_dc[t] = 0.0;
    } else if (want_grad) {
      double da, db; pcbvn::pbvn_grad_ab(a, b, rho, da, db);
      const double dr = pcbvn::dbvn(a, b, rho);
      dlogP_dd[0] = da / (s1 * P);
      dlogP_dd[1] = db / (s2 * P);
      for (int t = 0; t < nc; ++t) {
        const double ds1 = dOm[t][0] / (2.0 * s1), ds2 = dOm[t][3] / (2.0 * s2);
        const double drho = dOm[t][1] / (s1 * s2) - rho * (ds1 / s1 + ds2 / s2);
        const double dA = -a * ds1 / s1, dB = -b * ds2 / s2;
        dlogP_dc[t] = (da * dA + db * dB + dr * drho) / P;
      }
    }
    return lp;
  }
  // GHK
  std::vector<double> C;
  if (!chol_small(mc, Om, C)) return R_NegInf;      // a singular covariance: the case has no finite log-probability
  const int nd = want_grad ? mc + nc : 0;          // derivative directions: d_t, then c_t
  std::vector< std::vector<double> > dC(nc);
  if (want_grad) for (int t = 0; t < nc; ++t) chol_diff(mc, C, dOm[t], dC[t]);
  const int R = s.R;
  std::vector<double> b(mc), e(mc), db(nd * mc), de(nd * mc), dlogPr(nd), gsum(nd, 0.0);
  double Psum = 0.0;
  for (int r = 0; r < R; ++r) {
    const long long idx = 100 + (long long) i * R + r + 1;
    double logPr = 0.0;
    if (want_grad) std::fill(dlogPr.begin(), dlogPr.end(), 0.0);
    for (int st = 0; st < mc; ++st) {
      double num = d[st];
      for (int l = 0; l < st; ++l) num -= C[st * mc + l] * e[l];
      b[st] = num / C[st * mc + st];
      const double lPhi = R::pnorm(b[st], 0.0, 1.0, 1, 1);
      logPr += lPhi;
      double lam = 0.0, ratio = 0.0, u = 0.0;
      if (want_grad) lam = std::exp(R::dnorm(b[st], 0.0, 1.0, 1) - lPhi);
      if (st < mc - 1) {
        u = radical_inverse(idx, PRIMES[st]);
        double target = u * std::exp(lPhi);
        if (target < 1e-300) target = 1e-300;
        e[st] = R::qnorm(target, 0.0, 1.0, 1, 0);
        if (want_grad) ratio = std::exp(R::dnorm(b[st], 0.0, 1.0, 1) - R::dnorm(e[st], 0.0, 1.0, 1));
      }
      if (want_grad) {
        for (int v = 0; v < nd; ++v) {
          const bool is_d = v < mc;
          const std::vector<double>* dCv = is_d ? nullptr : &dC[v - mc];
          double val = is_d ? ((v == st) ? 1.0 : 0.0) : 0.0;
          for (int l = 0; l < st; ++l) {
            val -= C[st * mc + l] * de[v * mc + l];
            if (!is_d) val -= (*dCv)[st * mc + l] * e[l];
          }
          val /= C[st * mc + st];
          if (!is_d) val -= b[st] * (*dCv)[st * mc + st] / C[st * mc + st];
          db[v * mc + st] = val;
          dlogPr[v] += lam * val;
          if (st < mc - 1) de[v * mc + st] = u * ratio * val;
        }
      }
    }
    const double Pr = std::exp(logPr);
    Psum += Pr;
    if (want_grad) for (int v = 0; v < nd; ++v) gsum[v] += Pr * dlogPr[v];
  }
  const double P = std::max(Psum / R, 1e-300);
  if (want_grad) {
    const bool floored = !(Psum / R > 1e-300);                       // flat at the floor, as for J = 3
    for (int t = 0; t < mc; ++t) dlogP_dd[t] = floored ? 0.0 : (gsum[t] / R) / P;
    for (int t = 0; t < nc; ++t) dlogP_dc[t] = floored ? 0.0 : (gsum[mc + t] / R) / P;
  }
  return std::log(P);
}

void mnp_case(int i, const PSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
              const NumericMatrix& off, const IntegerMatrix& avail, int yi, bool want_grad, bool all_probs, PCase& c) {
  const int J = s.J;
  c.V.assign(J, R_NegInf); c.P.assign(J, 0.0); c.dV.assign(J, 0.0); c.dc.assign(s.nc, 0.0);
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    double v = off(i, j);
    if (s.bpos[j] >= 0) for (int cc = 0; cc < s.p; ++cc) v += X(i, cc) * theta[s.bpos[j] + cc];
    for (int cc = 0; cc < s.q; ++cc) v += W(i * J + j, cc) * theta[s.gpos + cc];
    c.V[j] = v;
  }
  std::vector<double> L; std::vector< std::vector<double> > dL;
  build_L(s, theta.begin() + s.cpos, L, dL);                 // one past the end when there are no covariance parameters (J = 2): not dereferenced
  std::vector<int> comp; std::vector<double> dd, dc;
  if (all_probs) for (int j = 0; j < J; ++j) if (avail(i, j)) c.P[j] = std::exp(mnp_logprob(i, j, s, c.V, L, dL, avail, false, comp, dd, dc));
  c.lli = mnp_logprob(i, yi, s, c.V, L, dL, avail, want_grad, comp, dd, dc);
  c.P[yi] = std::exp(c.lli);
  if (!want_grad) return;
  for (size_t t = 0; t < comp.size(); ++t) { c.dV[yi] += dd[t]; c.dV[comp[t]] -= dd[t]; }
  for (int t = 0; t < s.nc; ++t) c.dc[t] = dc[t];
}

void mnp_score(int i, const PSpec& s, const PCase& c, const NumericMatrix& X, const NumericMatrix& W,
               const IntegerMatrix& avail, std::vector<double>& out) {
  std::fill(out.begin(), out.end(), 0.0);
  const int J = s.J;
  for (int k = 0; k < J; ++k) {
    if (!avail(i, k)) continue;
    const double dV = c.dV[k];
    if (s.bpos[k] >= 0) for (int cc = 0; cc < s.p; ++cc) out[s.bpos[k] + cc] += dV * X(i, cc);
    for (int cc = 0; cc < s.q; ++cc) out[s.gpos + cc] += dV * W(i * J + k, cc);
  }
  for (int t = 0; t < s.nc; ++t) out[s.cpos + t] = c.dc[t];
}

} // namespace

// The loops over cases below check every 100 cases for a user interrupt (and
// so for R's time limits): one evaluation over a large sample takes about a
// second.
// [[Rcpp::export]]
NumericVector pc_mnp_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                             IntegerVector y, int J, int base, int R) {
  PSpec s = mnp_spec(theta, X, W, y, J, base, R);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  PCase c; NumericVector out(s.n);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    mnp_case(i, s, theta, X, W, off, avail, y[i], false, false, c); out[i] = c.lli;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pc_mnp_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                              IntegerVector y, NumericVector w, int J, int base, int R) {
  PSpec s = mnp_spec(theta, X, W, y, J, base, R);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  if (w.size() != s.n) stop("internal: weights do not match the cases.");
  PCase c; std::vector<double> sc(s.K); NumericVector g(s.K);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    mnp_case(i, s, theta, X, W, off, avail, y[i], true, false, c);
    mnp_score(i, s, c, X, W, avail, sc);
    for (int k = 0; k < s.K; ++k) g[k] += w[i] * sc[k];
  }
  return g;
}

// [[Rcpp::export]]
NumericMatrix pc_mnp_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                                IntegerVector y, int J, int base, int R) {
  PSpec s = mnp_spec(theta, X, W, y, J, base, R);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  PCase c; std::vector<double> sc(s.K); NumericMatrix S(s.n, s.K);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    mnp_case(i, s, theta, X, W, off, avail, y[i], true, false, c);
    mnp_score(i, s, c, X, W, avail, sc);
    for (int k = 0; k < s.K; ++k) S(i, k) = sc[k];
  }
  return S;
}

// probabilities of every alternative (GHK probabilities need not sum exactly to one) and the utilities
// [[Rcpp::export]]
List pc_mnp_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                      IntegerVector y, int J, int base, int R) {
  PSpec s = mnp_spec(theta, X, W, y, J, base, R);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  PCase c; NumericMatrix P(s.n, J), Vm(s.n, J);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    int yi = y[i];
    if (yi < 0 || yi >= J || !avail(i, yi)) { yi = 0; while (yi < J && !avail(i, yi)) ++yi; if (yi == J) yi = 0; }
    mnp_case(i, s, theta, X, W, off, avail, yi, false, true, c);
    for (int j = 0; j < J; ++j) { P(i, j) = c.P[j]; Vm(i, j) = c.V[j]; }
  }
  return List::create(_["prob"] = P, _["V"] = Vm);
}

// the covariance of the differences from the base implied by the structural block (for reporting / tests)
// [[Rcpp::export]]
NumericMatrix pc_mnp_sigma_cpp(NumericVector cpar, int m) {
  std::vector<double> L(m * m, 0.0); L[0] = 1.0;
  int t = 0;
  for (int r = 0; r < m; ++r) for (int c = 0; c <= r; ++c) if (!(r == 0 && c == 0)) { L[r * m + c] = (r == c) ? std::exp(cpar[t]) : cpar[t]; ++t; }
  NumericMatrix S(m, m);
  for (int a = 0; a < m; ++a) for (int b = 0; b < m; ++b) { double v = 0.0; for (int c = 0; c < m; ++c) v += L[a * m + c] * L[b * m + c]; S(a, b) = v; }
  return S;
}
