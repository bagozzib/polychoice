// pc_bvn.h -- bivariate standard-normal probabilities and their partial
// derivatives, used by the three-alternative multinomial probit (mnp_ll.cpp).
// Ported from the author's iop package (src/iop_bvn.h, same algorithm and
// validation): the rectangle probability of Drezner & Wesolowsky (1990) as
// refined by Genz (2004), Gauss-Legendre quadrature of the integral
// representation in the correlation with 6/12/20 points depending on |rho|
// and a separate expansion for |rho| > 0.925; deterministic with absolute
// error around 5e-16.
//
// References: Drezner, Z. and Wesolowsky, G.O. (1990), "On the computation of
// the bivariate normal integral", J. Statist. Comput. Simul. 35, 101-107.
// Genz, A. (2004), "Numerical computation of rectangular bivariate and
// trivariate normal and t probabilities", Statistics and Computing 14, 251-260.
#ifndef PC_BVN_H
#define PC_BVN_H

#include <Rmath.h>
#include <cmath>
#include <algorithm>

namespace pcbvn {

inline double Phi(double x)  { return R::pnorm(x, 0.0, 1.0, 1, 0); }
inline double phi(double x)  { return R::dnorm(x, 0.0, 1.0, 0); }

static const double GL_X6[3]  = {0.2386191860831971, 0.6612093864662647, 0.9324695142031519};
static const double GL_W6[3]  = {0.4679139345726905, 0.3607615730481394, 0.1713244923791702};
static const double GL_X12[6] = {0.1252334085114698, 0.3678314989981801, 0.5873179542866176,
                                 0.7699026741943045, 0.9041172563704747, 0.9815606342467192};
static const double GL_W12[6] = {0.2491470458134025, 0.2334925365383561, 0.2031674267230653,
                                 0.1600783285433461, 0.1069393259953181, 0.0471753363865117};
static const double GL_X20[10]= {0.0765265211334978, 0.2277858511416454, 0.3737060887154197,
                                 0.5108670019508271, 0.6360536807265152, 0.7463319064601511,
                                 0.8391169718222189, 0.9122344282513261, 0.9639719272779138,
                                 0.9931285991850950};
static const double GL_W20[10]= {0.1527533871307258, 0.1491729864726030, 0.1420961093183802,
                                 0.1316886384491760, 0.1181945319615197, 0.1019301198172408,
                                 0.0832767415767051, 0.0626720483341093, 0.0406014298003871,
                                 0.0176140071391521};

// Upper-orthant probability P(X > dh, Y > dk), correlation r, finite dh, dk.
inline double bvnu(double dh, double dk, double r) {
  const double TWOPI = 6.283185307179586;
  const double *X; const double *W; int LG;
  if (std::fabs(r) < 0.3)       { X = GL_X6;  W = GL_W6;  LG = 3; }
  else if (std::fabs(r) < 0.75) { X = GL_X12; W = GL_W12; LG = 6; }
  else                          { X = GL_X20; W = GL_W20; LG = 10; }
  double h = dh, k = dk, hk = h * k, bvn = 0.0;
  if (std::fabs(r) < 0.925) {
    if (std::fabs(r) > 0.0) {
      double hs = (h * h + k * k) / 2.0, asr = std::asin(r);
      for (int i = 0; i < LG; i++) {
        for (int is = -1; is <= 1; is += 2) {
          double sn = std::sin(asr * (is * X[i] + 1.0) / 2.0);
          bvn += W[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
        }
      }
      bvn = bvn * asr / (2.0 * TWOPI);
    }
    bvn += Phi(-h) * Phi(-k);
  } else {
    if (r < 0.0) { k = -k; hk = -hk; }
    if (std::fabs(r) < 1.0) {
      double as = (1.0 - r) * (1.0 + r), a = std::sqrt(as), bs = (h - k) * (h - k);
      double c = (4.0 - hk) / 8.0, d = (12.0 - hk) / 16.0, asr = -(bs / as + hk) / 2.0;
      if (asr > -100.0)
        bvn = a * std::exp(asr) * (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);
      if (-hk < 100.0) {
        double b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2.0) * std::sqrt(TWOPI) * Phi(-b / a) * b * (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }
      a = a / 2.0;
      for (int i = 0; i < LG; i++) {
        for (int is = -1; is <= 1; is += 2) {
          double xs = a * (is * X[i] + 1.0); xs = xs * xs;
          double rs = std::sqrt(1.0 - xs), asr2 = -(bs / xs + hk) / 2.0;
          if (asr2 > -100.0)
            bvn += a * W[i] * std::exp(asr2) *
                   (std::exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs - (1.0 + c * xs * (1.0 + d * xs)));
        }
      }
      bvn = -bvn / TWOPI;
    }
    if (r > 0.0) {
      bvn += Phi(-std::max(h, k));
    } else {
      bvn = -bvn;
      if (k > h) bvn += Phi(k) - Phi(h);
    }
  }
  if (bvn < 0.0) bvn = 0.0;
  if (bvn > 1.0) bvn = 1.0;
  return bvn;
}

// Lower-tail rectangle probability Phi2(a, b; r) = P(X <= a, Y <= b).
inline double pbvn(double a, double b, double r) {
  if (ISNAN(a) || ISNAN(b) || ISNAN(r)) return NA_REAL;
  if (a == R_NegInf || b == R_NegInf) return 0.0;
  if (a == R_PosInf) return (b == R_PosInf) ? 1.0 : Phi(b);
  if (b == R_PosInf) return Phi(a);
  if (r >= 1.0)  return Phi(std::min(a, b));
  if (r <= -1.0) return std::max(0.0, Phi(a) + Phi(b) - 1.0);
  return bvnu(-a, -b, r);
}

// Bivariate standard-normal density at (a, b) with correlation r (|r| < 1):
// the partial derivative of Phi2(a, b; r) with respect to r.
inline double dbvn(double a, double b, double r) {
  double s2 = 1.0 - r * r;
  return std::exp(-(a * a - 2.0 * r * a * b + b * b) / (2.0 * s2)) / (6.283185307179586 * std::sqrt(s2));
}

// Partial derivatives of Phi2(a, b; r) with respect to a and b.
inline void pbvn_grad_ab(double a, double b, double r, double &da, double &db) {
  da = db = 0.0;
  if (a == R_NegInf || b == R_NegInf) return;
  if (a == R_PosInf && b == R_PosInf) return;
  if (b == R_PosInf) { da = phi(a); return; }
  if (a == R_PosInf) { db = phi(b); return; }
  double s = std::sqrt(std::max(1.0 - r * r, 1e-300));
  da = phi(a) * Phi((b - r * a) / s);
  db = phi(b) * Phi((a - r * b) / s);
}

} // namespace pcbvn
#endif
