// Bottom layer of the summary-statistics PRS methods (viprs.cpp, sbayesrc.cpp): one stored LD row and the two
// operations every coordinate update needs on it, for float32 and IEEE float16 storage.
//   put(row, i, x): store x (rounded to the storage type), return the stored value as float
//   get(row, i):    stored value as float
//   dot(row, s, k): sum_i row_i s_i (float lanes, double total)
//   axpy(s, row, a, k): s += a row
// Ported from viprs_eigen/src/_engine.cpp (Store<T>, h2f, f2h). dot and axpy run on the SIMD kernels of the package
// dispatch (la_dispatch.h: AVX2 + FMA + F16C, the original's build, or the generic copy).
#ifndef PLINK4R_VI_STORE_H
#define PLINK4R_VI_STORE_H

#include <cstdint>
#include "half.h"
#include "la_dispatch.h"

namespace plink4r {

template <typename T> struct RowStore;

template <> struct RowStore<float> {
  static inline float put(float* row, int i, float x) { row[i] = x; return x; }
  static inline float get(const float* row, int i) { return row[i]; }
  static inline float dot(const float* v, const float* s, int k) { return la().vi_dot_f32(v, s, k); }
  static inline void axpy(float* s, const float* v, float a, int k) { la().vi_axpy_f32(s, v, a, k); }
};

template <> struct RowStore<uint16_t> {
  static inline float put(uint16_t* row, int i, float x) {
    const uint16_t h = float_to_half(x);
    row[i] = h;
    return half_to_float(h);
  }
  static inline float get(const uint16_t* row, int i) { return half_to_float(row[i]); }
  static inline float dot(const uint16_t* v, const float* s, int k) { return la().vi_dot_f16(v, s, k); }
  static inline void axpy(float* s, const uint16_t* v, float a, int k) { la().vi_axpy_f16(s, v, a, k); }
};

}  // namespace plink4r
#endif
