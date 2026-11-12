// hev_ll.cpp -- the heteroskedastic extreme value model (Bhat 1995): independent
// Gumbel errors with alternative-specific scales theta_j (the base alternative's
// scale is 1), so that IIA is relaxed through the scales.
//
//   V_ij as in mnl_ll.cpp. With w = eps_j / theta_j standard Gumbel,
//     P(j) = int_{-inf}^{inf} exp(-S_j(w)) f(w) dw,   f(w) = exp(-w - e^{-w}),
//     S_j(w) = sum_{k != j, available} h_k(w),  h_k(w) = exp(-a_k - r_k w),
//     a_k = (V_j - V_k)/theta_k,  r_k = theta_j/theta_k.
//   With every theta = 1, S_j(w) = (rho - 1) e^{-w} with rho = 1 + sum_k exp(-a_k)
//   and P(j) = 1/rho (the multinomial logit).
//
// Quadrature. For each case and alternative the integrand g(w) = exp(-S(w)) f(w)
// is log-concave (d^2 log g / dw^2 = -sum_k r_k^2 h_k - e^{-w} < 0), so its mode
// w* > 0 is found by bisection on d log g / dw = sum_k r_k h_k - 1 + e^{-w}, and
// the window is where log g has fallen 36 below its maximum (a relative 1e-16;
// by concavity nothing lies beyond), at most 200 units of w to either side.
// Two rules share that window:
// * r_max <= 20 (every competitor's scale at least a twentieth of theta_j): the
//   integrand is analytic with features no narrower than 1/r_max, so the
//   trapezoidal rule with step 0.35 / max(1, r_max) converges spectrally (error
//   about exp(-pi^2 / (r_max h))); at most 400 / h nodes.
// * r_max > 20: exp(-h_k(w)) for a competitor with a much smaller scale is a
//   near step of width 1/r_k at w = -a_k/r_k, which a uniform rule needs O(r_k)
//   nodes to resolve. The window is then covered by Gauss-Kronrod (7-15) panels
//   of at most 5 units, bisecting the panel with the largest error estimate
//   until the summed error of g / g(w*) and of the two companions that carry
//   the derivatives, g sum_k r_k h_k and g sum_k r_k h_k |a_k + r_k w|, is
//   below 1e-10 of their values (at most 2048 panels, 50 levels): O(log r_k)
//   panels per step. The initial panels carry breakpoints around each sharp
//   competitor's step (see below), and the error companion of the scale
//   derivative uses sqrt(1 + e^2) rather than |e|, whose kink sat at the peak
//   of the derivative spike and slowed the refinement there. Checked against
//   R's integrate() to the test tolerance of 1e-9 for scale ratios up to
//   1.6e5 (test-families.R).
// The first implementation used the trapezoidal rule alone with the ratio in
// the step capped at 60: capped, it under-resolved the step and biased the
// likelihood upward on TravelMode; uncapped, an optimizer excursion asked for
// 1e11 nodes. log P is returned, so the log-likelihood stays finite where P
// underflows. Bhat's Gauss-Laguerre scheme in u = e^{-w} was tried first and
// converges only algebraically when the scale ratios differ from 1 (u^{r} has a
// singular derivative at 0 for r < 1), which is why it is not used.
//
// Parameter vector theta: beta (p (J-1)), gamma (q), then t_j = log theta_j for
// the non-base alternatives in level order.
//
// Derivatives of l = log P(y), with pi_q = omega_q g(w_q) / P the normalized
// quadrature weights:
//   d l / d V_y     =  sum_q pi_q sum_k h_kq / theta_k
//   d l / d V_k     = -sum_q pi_q h_kq / theta_k                                   (k != y)
//   d l / d theta_y =  sum_q pi_q w_q sum_k h_kq / theta_k
//   d l / d theta_k = -sum_q pi_q (h_kq / theta_k) (a_k + r_k w_q)                   (k != y)
//   d l / d t = theta d l / d theta. The nodes and weights are held fixed in the
//   derivatives (the integral does not depend on them).
// Validated against numDeriv in tests/testthat/test-families.R; the Hessian is
// taken by finite differences of this gradient (R side).
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>
using namespace Rcpp;

namespace {

struct HSpec {
  int n, J, p, q, base, K;
  std::vector<int> bpos, tpos;   // tpos[j]: index of t_j in theta, -1 for the base
  int gpos;
};

HSpec hev_spec(const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W, const IntegerVector& y, int J, int base) {
  HSpec s;
  s.n = y.size(); s.J = J; s.p = X.ncol(); s.q = W.ncol(); s.base = base;
  s.bpos.assign(J, -1); s.tpos.assign(J, -1);
  int pos = 0;
  for (int j = 0; j < J; ++j) if (j != base) { s.bpos[j] = pos; pos += s.p; }
  s.gpos = pos; pos += s.q;
  for (int j = 0; j < J; ++j) if (j != base) { s.tpos[j] = pos; pos += 1; }
  s.K = pos;
  if (theta.size() != s.K) stop("internal: theta has length %d, expected %d.", (int) theta.size(), s.K);
  if (X.nrow() != s.n || W.nrow() != s.n * J) stop("internal: design rows do not match the response.");
  return s;
}

struct HCase {
  std::vector<double> V, th, P;    // utilities, scales, probabilities
  std::vector<double> dV, dth;     // d l / d V_k, d l / d theta_k for the chosen alternative
  double lli;
};

// Gauss-Kronrod 7-15 (QUADPACK qk15): abscissae on [0, 1] in decreasing order,
// the Kronrod weights, and the Gauss weights of the odd abscissae and the centre
const double GK_X[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                        0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                        0.586087235467691130294144838258730, 0.405845151377397166906606412076961,
                        0.207784955007898467600689403773245, 0.0};
const double GK_WK[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                         0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                         0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                         0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
const double GK_WG[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                         0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

// S(w) and sum_k r_k h_k(w) = -S'(w) over the competitors of one alternative
inline void hev_sums(const std::vector<double>& a, const std::vector<double>& r, double w, double& S, double& rS) {
  S = 0.0; rS = 0.0;
  for (size_t t = 0; t < a.size(); ++t) {
    const double e = -a[t] - r[t] * w;
    const double hv = (e > 700.0) ? INFINITY : std::exp(e);
    S += hv; rS += r[t] * hv;
  }
}
inline double hev_logg(const std::vector<double>& a, const std::vector<double>& r, double w) {
  double S, rS; hev_sums(a, r, w, S, rS);
  return -S - w - std::exp(-w);
}
inline double hev_dlogg(const std::vector<double>& a, const std::vector<double>& r, double w) {
  double S, rS; hev_sums(a, r, w, S, rS);
  return rS - 1.0 + std::exp(-w);
}

// the quadrature nodes of one (case, alternative): positions, weights, the
// integrand scaled by its maximum, and (when asked) the h_k values node-major
struct HNodes {
  std::vector<double> w, wt, f, hk;
  void add(double w_, double wt_, double f_, const double* h, int m) {
    w.push_back(w_); wt.push_back(wt_); f.push_back(f_);
    if (h) for (int k = 0; k < m; ++k) hk.push_back(h[k]);
  }
};

// the scaled integrand exp(log g(w) - Ls) and the derivative companions at w
inline double hev_G(const std::vector<double>& a, const std::vector<double>& r, double w, double Ls, double* h, double& A, double& B) {
  double S = 0.0; A = 0.0; B = 0.0;
  for (size_t k = 0; k < a.size(); ++k) {
    const double e = -a[k] - r[k] * w;
    const double hv = (e > 700.0) ? INFINITY : std::exp(e);
    if (h) h[k] = hv;
    S += hv; A += r[k] * hv; B += r[k] * hv * std::sqrt(1.0 + e * e);
  }
  const double lg = -S - w - std::exp(-w);
  const double f = std::isfinite(lg) ? std::exp(lg - Ls) : 0.0;
  if (f == 0.0) { A = 0.0; B = 0.0; } else { A *= f; B *= f; }
  return f;
}

// integrals whose adaptive refinement stopped at the panel or depth limit
// before reaching its tolerance (read and reset by pc_hev_unsettled)
static int hev_unsettled_count = 0;

struct HPanel {
  double a, b; int depth;
  double x[15], f[15];
  std::vector<double> hk;            // 15 * m when kept
  double k15[3], err[3];
};

void hev_panel(HPanel& p, const std::vector<double>& a, const std::vector<double>& r, double Ls, bool keep_h) {
  const int m = a.size();
  const double cc = 0.5 * (p.a + p.b), hl = 0.5 * (p.b - p.a);
  if (keep_h) p.hk.assign(15 * m, 0.0);
  for (int t = 0; t < 7; ++t) { p.x[t] = cc - hl * GK_X[t]; p.x[14 - t] = cc + hl * GK_X[t]; }
  p.x[7] = cc;
  double k15[3] = {0.0, 0.0, 0.0}, g7[3] = {0.0, 0.0, 0.0};
  for (int t = 0; t < 15; ++t) {
    double A, B;
    p.f[t] = hev_G(a, r, p.x[t], Ls, keep_h ? &p.hk[t * m] : 0, A, B);
    const double v[3] = {p.f[t], A, B};
    const int u = (t <= 7) ? t : 14 - t;
    for (int c = 0; c < 3; ++c) {
      k15[c] += GK_WK[u] * v[c];
      if (u == 1 || u == 3 || u == 5) g7[c] += GK_WG[(u - 1) / 2] * v[c];
      else if (u == 7) g7[c] += GK_WG[3] * v[c];
    }
  }
  for (int c = 0; c < 3; ++c) { p.k15[c] = hl * k15[c]; p.err[c] = std::fabs(hl * (k15[c] - g7[c])); }
}

// log P(j) and (when want_grad) the derivatives of log P(j) for one case and alternative j
double hev_logprob(int i, int j, const HSpec& s, const HCase& c, const IntegerMatrix& avail,
                   bool want_grad, std::vector<double>& dV, std::vector<double>& dth) {
  const int J = s.J;
  std::vector<int> others;
  for (int k = 0; k < J; ++k) if (k != j && avail(i, k)) others.push_back(k);
  const int m = others.size();
  if (want_grad) { std::fill(dV.begin(), dV.end(), 0.0); std::fill(dth.begin(), dth.end(), 0.0); }
  if (m == 0) return 0.0;                                  // the only available alternative: P = 1
  std::vector<double> a(m), r(m);
  double rmax = 0.0;
  for (int t = 0; t < m; ++t) {
    const int k = others[t];
    a[t] = (c.V[j] - c.V[k]) / c.th[k]; r[t] = c.th[j] / c.th[k];
    rmax = std::max(rmax, r[t]);
  }
  // the mode of the log-integrand: L'(0) = sum_k r_k h_k(0) > 0; double the
  // upper bracket until L' < 0, then bisect
  double wa = 0.0, wb = 1.0;
  while (hev_dlogg(a, r, wb) > 0.0) { wb *= 2.0; if (wb > 1e15) return -INFINITY; }
  for (int it = 0; it < 40; ++it) { const double wm = 0.5 * (wa + wb); if (hev_dlogg(a, r, wm) > 0.0) wa = wm; else wb = wm; }
  const double ws = 0.5 * (wa + wb), Ls = hev_logg(a, r, ws);
  if (!std::isfinite(Ls)) return -INFINITY;
  const double drop = 36.0, wmax = 200.0;
  HNodes nd;
  std::vector<double> hbuf(m);
  double total = 0.0;
  if (rmax <= 20.0) {
    // trapezoidal rule from the mode outwards until the log-integrand has
    // dropped by 36 (or 200 units); nodes are filled during the scan
    const double h = 0.35 / std::max(1.0, rmax);
    double A, B;
    double f0 = hev_G(a, r, ws, Ls, want_grad ? &hbuf[0] : 0, A, B);
    nd.add(ws, h, f0, want_grad ? &hbuf[0] : 0, m); total += h * f0;
    for (int side = -1; side <= 1; side += 2) {
      for (int q = 1; q * h <= wmax; ++q) {
        const double w = ws + side * q * h;
        const double f = hev_G(a, r, w, Ls, want_grad ? &hbuf[0] : 0, A, B);
        nd.add(w, h, f, want_grad ? &hbuf[0] : 0, m); total += h * f;
        if (f < std::exp(-drop)) break;
      }
    }
  } else {
    // adaptive Gauss-Kronrod on the window
    double Lw = ws, Uw = ws;
    for (int k = 1; k <= 200; ++k) { Lw = ws - k; if (hev_logg(a, r, Lw) < Ls - drop) break; }
    for (int k = 1; k <= 200; ++k) { Uw = ws + k; if (hev_logg(a, r, Uw) < Ls - drop) break; }
    const double tol = 1e-10;
    const int max_panels = 2048, max_depth = 50;
    // initial panels: the window in pieces of at most 5 units plus breakpoints
    // around the step of every competitor with r_k > 1, where exp(-h_k) rises
    // from 0 to 1 over w*_k - 4 / r_k .. w*_k + 40 / r_k, w*_k = -a_k / r_k
    // (h_k = exp(-r_k (w - w*_k)) runs from e^4 to e^-40 there). The breakpoints
    // move continuously with theta and resolve the transition at its own scale.
    // Without them a bisection could leave the transition in a sliver at the
    // end of a panel where no node falls: both rules missed it, the error
    // estimate was small, the integral came out short by up to 1e-5, and the
    // log-likelihood and its gradient jumped as theta moved, which corrupted the
    // finite-difference Hessian (TravelMode boundary fit, ledger 2026-10-14).
    std::vector<double> bp;
    const int n0 = std::max(1, (int) std::ceil((Uw - Lw) / 5.0));
    for (int q = 0; q <= n0; ++q) bp.push_back(Lw + (Uw - Lw) * q / n0);
    const double offs[7] = {-4.0, -1.0, 0.0, 1.0, 4.0, 12.0, 40.0};
    for (int t = 0; t < m; ++t) if (r[t] > 1.0) {
      const double wst = -a[t] / r[t];
      for (int o = 0; o < 7; ++o) { const double b = wst + offs[o] / r[t]; if (b > Lw && b < Uw) bp.push_back(b); }
    }
    std::sort(bp.begin(), bp.end());
    std::vector<double> cut;
    for (size_t q = 0; q < bp.size(); ++q)
      if (cut.empty() || bp[q] - cut.back() > 1e-12 * std::max(1.0, std::fabs(bp[q]))) cut.push_back(bp[q]);
    if (cut.size() < 2) cut.push_back(Uw); else cut.back() = Uw;
    std::vector<HPanel> panels;
    for (size_t q = 0; q + 1 < cut.size(); ++q) {
      HPanel p; p.a = cut[q]; p.b = cut[q + 1]; p.depth = 0;
      hev_panel(p, a, r, Ls, want_grad); panels.push_back(p);
    }
    for (;;) {
      double tot[3] = {0.0, 0.0, 0.0}, se[3] = {0.0, 0.0, 0.0};
      for (size_t q = 0; q < panels.size(); ++q) for (int cc = 0; cc < 3; ++cc) { tot[cc] += panels[q].k15[cc]; se[cc] += panels[q].err[cc]; }
      bool settled = true;
      for (int cc = 0; cc < 3; ++cc) if (se[cc] > tol * tot[cc]) settled = false;
      if (settled) break;
      if ((int) panels.size() >= max_panels) { if (hev_unsettled_count < INT_MAX) ++hev_unsettled_count; break; }
      size_t worst = 0; double worst_e = -1.0;
      for (size_t q = 0; q < panels.size(); ++q) {
        double e = 0.0;
        for (int cc = 0; cc < 3; ++cc) if (tot[cc] > 0.0) e += panels[q].err[cc] / tot[cc];
        if (e > worst_e) { worst_e = e; worst = q; }
      }
      if (panels[worst].depth >= max_depth) { if (hev_unsettled_count < INT_MAX) ++hev_unsettled_count; break; }
      HPanel left, right;
      left.a = panels[worst].a; right.b = panels[worst].b;
      left.b = right.a = 0.5 * (left.a + right.b);
      left.depth = right.depth = panels[worst].depth + 1;
      hev_panel(left, a, r, Ls, want_grad); hev_panel(right, a, r, Ls, want_grad);
      panels[worst] = left; panels.push_back(right);
    }
    for (size_t q = 0; q < panels.size(); ++q) {
      const HPanel& p = panels[q];
      const double hl = 0.5 * (p.b - p.a);
      total += p.k15[0];
      for (int t = 0; t < 15; ++t) { const int u = (t <= 7) ? t : 14 - t; nd.add(p.x[t], hl * GK_WK[u], p.f[t], want_grad ? &p.hk[t * m] : 0, m); }
    }
  }
  if (!(total > 0.0)) return -INFINITY;
  const double logP = Ls + std::log(total);
  if (!want_grad) return logP;
  for (size_t q = 0; q < nd.w.size(); ++q) {
    const double pi = nd.wt[q] * nd.f[q] / total;
    if (pi == 0.0) continue;
    const double w = nd.w[q];
    double sumh = 0.0;
    for (int k = 0; k < m; ++k) {
      const double hv = nd.hk[q * m + k];
      if (!std::isfinite(hv)) continue;
      const int kk = others[k];
      const double hk_th = hv / c.th[kk];
      sumh += hk_th;
      dV[kk] -= pi * hk_th;
      dth[kk] -= pi * hk_th * (a[k] + r[k] * w);
    }
    dV[j] += pi * sumh;
    dth[j] += pi * w * sumh;
  }
  return logP;
}

void hev_case(int i, const HSpec& s, const NumericVector& theta, const NumericMatrix& X, const NumericMatrix& W,
              const NumericMatrix& off, const IntegerMatrix& avail, int yi, bool want_grad, bool all_probs, HCase& c) {
  const int J = s.J;
  c.V.assign(J, R_NegInf); c.th.assign(J, 1.0); c.P.assign(J, 0.0);
  c.dV.assign(J, 0.0); c.dth.assign(J, 0.0);
  for (int j = 0; j < J; ++j) if (s.tpos[j] >= 0) c.th[j] = std::exp(theta[s.tpos[j]]);
  for (int j = 0; j < J; ++j) {
    if (!avail(i, j)) continue;
    double v = off(i, j);
    if (s.bpos[j] >= 0) for (int cc = 0; cc < s.p; ++cc) v += X(i, cc) * theta[s.bpos[j] + cc];
    for (int cc = 0; cc < s.q; ++cc) v += W(i * J + j, cc) * theta[s.gpos + cc];
    c.V[j] = v;
  }
  std::vector<double> dV(J), dth(J);
  if (all_probs) for (int j = 0; j < J; ++j) if (avail(i, j)) c.P[j] = std::exp(hev_logprob(i, j, s, c, avail, false, dV, dth));
  const double lPy = hev_logprob(i, yi, s, c, avail, want_grad, c.dV, c.dth);
  c.P[yi] = std::exp(lPy);
  c.lli = lPy;
}

void hev_score(int i, const HSpec& s, const HCase& c, const NumericMatrix& X, const NumericMatrix& W,
               const IntegerMatrix& avail, std::vector<double>& out) {
  std::fill(out.begin(), out.end(), 0.0);
  const int J = s.J;
  for (int k = 0; k < J; ++k) {
    if (!avail(i, k)) continue;
    const double dV = c.dV[k];
    if (s.bpos[k] >= 0) for (int cc = 0; cc < s.p; ++cc) out[s.bpos[k] + cc] += dV * X(i, cc);
    for (int cc = 0; cc < s.q; ++cc) out[s.gpos + cc] += dV * W(i * J + k, cc);
    if (s.tpos[k] >= 0) out[s.tpos[k]] += c.dth[k] * c.th[k];       // d/dt = theta d/dtheta
  }
}

} // namespace

// the number of HEV integrals that stopped short of their tolerance since the
// last reset; `reset` sets the counter to 0 after reading it
// [[Rcpp::export]]
int pc_hev_unsettled(bool reset) {
  const int k = hev_unsettled_count;
  if (reset) hev_unsettled_count = 0;
  return k;
}

// The loops over cases below check every 100 cases for a user interrupt (and
// so for R's time limits): one evaluation over a large sample takes about a
// second, and an integral at its panel or depth limit takes longer.
// [[Rcpp::export]]
NumericVector pc_hev_lli_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                             IntegerVector y, int J, int base) {
  HSpec s = hev_spec(theta, X, W, y, J, base);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  HCase c; NumericVector out(s.n);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    hev_case(i, s, theta, X, W, off, avail, y[i], false, false, c); out[i] = c.lli;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pc_hev_grad_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                              IntegerVector y, NumericVector w, int J, int base) {
  HSpec s = hev_spec(theta, X, W, y, J, base);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  if (w.size() != s.n) stop("internal: weights do not match the cases.");
  HCase c; std::vector<double> sc(s.K); NumericVector g(s.K);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    hev_case(i, s, theta, X, W, off, avail, y[i], true, false, c);
    hev_score(i, s, c, X, W, avail, sc);
    for (int k = 0; k < s.K; ++k) g[k] += w[i] * sc[k];
  }
  return g;
}

// [[Rcpp::export]]
NumericMatrix pc_hev_scores_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                                IntegerVector y, int J, int base) {
  HSpec s = hev_spec(theta, X, W, y, J, base);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  HCase c; std::vector<double> sc(s.K); NumericMatrix S(s.n, s.K);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    hev_case(i, s, theta, X, W, off, avail, y[i], true, false, c);
    hev_score(i, s, c, X, W, avail, sc);
    for (int k = 0; k < s.K; ++k) S(i, k) = sc[k];
  }
  return S;
}

// [[Rcpp::export]]
List pc_hev_probs_cpp(NumericVector theta, NumericMatrix X, NumericMatrix W, NumericMatrix off, IntegerMatrix avail,
                      IntegerVector y, int J, int base) {
  HSpec s = hev_spec(theta, X, W, y, J, base);
  if (off.nrow() != s.n || off.ncol() != J || avail.nrow() != s.n || avail.ncol() != J) stop("internal: offset or availability dimensions do not match.");
  HCase c; NumericMatrix P(s.n, J), Vm(s.n, J);
  for (int i = 0; i < s.n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    int yi = y[i];
    if (yi < 0 || yi >= J || !avail(i, yi)) { yi = 0; while (yi < J && !avail(i, yi)) ++yi; if (yi == J) yi = 0; }
    hev_case(i, s, theta, X, W, off, avail, yi, false, true, c);
    for (int j = 0; j < J; ++j) { P(i, j) = c.P[j]; Vm(i, j) = c.V[j]; }
  }
  return List::create(_["prob"] = P, _["V"] = Vm);
}
