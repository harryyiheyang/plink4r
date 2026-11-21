// Brent's one-dimensional minimiser on [a, b] (golden section with parabolic steps; Brent 1973, "localmin"),
// without derivatives. The package's one scalar optimiser.
#ifndef PLINK4R_BRENT_H
#define PLINK4R_BRENT_H

#include <cmath>
#include <limits>

namespace plink4r {

struct BrentResult {
  double x = 0, fx = 0;
  int evals = 0;
};

// Minimises f on [a, b]; stops when the bracket is within 2 (tol |x| + tol_abs) of x, when the three best points
// are distinct and their values within f_tol of each other (f_tol > 0), or after max_eval calls.
template <class F>
inline BrentResult brent_min(F f, double a, double b, double tol, double tol_abs, int max_eval, double f_tol = 0) {
  const double gold = 0.3819660112501051;  // (3 - sqrt(5)) / 2
  double x = a + gold * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0, e = 0;
  BrentResult out;
  out.evals = 1;
  while (out.evals < max_eval) {
    const double m = 0.5 * (a + b);
    const double t1 = tol * std::fabs(x) + tol_abs, t2 = 2 * t1;
    if (std::fabs(x - m) <= t2 - 0.5 * (b - a)) break;
    if (f_tol > 0 && v != w && w != x && v != x && fv - fx <= f_tol) break;
    bool golden = true;
    if (std::fabs(e) > t1) {  // parabola through (v, fv), (w, fw), (x, fx)
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p; else q = -q;
      r = e;
      e = d;
      if (std::fabs(p) < std::fabs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        const double u = x + d;
        if (u - a < t2 || b - u < t2) d = x < m ? t1 : -t1;
        golden = false;
      }
    }
    if (golden) {
      e = (x < m ? b : a) - x;
      d = gold * e;
    }
    const double u = std::fabs(d) >= t1 ? x + d : x + (d > 0 ? t1 : -t1);
    const double fu = f(u);
    ++out.evals;
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  out.x = x;
  out.fx = fx;
  return out;
}

}  // namespace plink4r
#endif
