// Bottom layer of the summary-statistics PRS methods (viprs.cpp, sbayesrc.cpp): one stored LD row and the two
// operations every coordinate update needs on it, for float32 and IEEE float16 storage.
//   put(row, i, x): store x (rounded to the storage type), return the stored value as float
//   get(row, i):    stored value as float
//   dot(row, s, k): sum_i row_i s_i (float lanes, double total)
//   axpy(s, row, a, k): s += a row
// Ported from viprs_eigen/src/_engine.cpp (Store<T>, h2f, f2h); F16C is used when the compiler targets it,
// otherwise the exact scalar binary16 conversion (round to nearest even).
#ifndef PLINK4R_VI_STORE_H
#define PLINK4R_VI_STORE_H

#include <cmath>
#include <cstdint>
#include <cstring>

#if defined(__F16C__)
#include <immintrin.h>
#endif

namespace plink4r {

inline float half_to_float(uint16_t h) {
#if defined(__F16C__)
  return _cvtsh_ss(h);
#else
  uint32_t sign = (uint32_t)(h & 0x8000) << 16, ex = (h >> 10) & 0x1f, man = h & 0x3ff, f;
  if (ex == 0) {
    if (man == 0) f = sign;
    else {
      ex = 127 - 15 + 1;
      while (!(man & 0x400)) { man <<= 1; --ex; }
      man &= 0x3ff;
      f = sign | (ex << 23) | (man << 13);
    }
  } else if (ex == 31) f = sign | 0x7f800000 | (man << 13);
  else f = sign | ((ex + 127 - 15) << 23) | (man << 13);
  float out;
  std::memcpy(&out, &f, 4);
  return out;
#endif
}

inline uint16_t float_to_half(float x) {
#if defined(__F16C__)
  return _cvtss_sh(x, _MM_FROUND_TO_NEAREST_INT);
#else
  uint32_t f;
  std::memcpy(&f, &x, 4);
  uint32_t sign = (f >> 16) & 0x8000;
  int32_t ex = (int32_t)((f >> 23) & 0xff) - 127 + 15;
  uint32_t man = f & 0x7fffff;
  if (((f >> 23) & 0xff) == 0xff) return (uint16_t)(sign | 0x7c00 | (man ? 0x200 : 0));
  if (ex >= 31) return (uint16_t)(sign | 0x7c00);
  if (ex <= 0) {
    if (ex < -10) return (uint16_t)sign;
    man |= 0x800000;
    int shift = 14 - ex;
    uint32_t half = man >> shift, rem = man & ((1u << shift) - 1), mid = 1u << (shift - 1);
    if (rem > mid || (rem == mid && (half & 1))) ++half;
    return (uint16_t)(sign | half);
  }
  uint32_t half = sign | ((uint32_t)ex << 10) | (man >> 13), rem = man & 0x1fff;
  if (rem > 0x1000 || (rem == 0x1000 && (half & 1))) ++half;
  return (uint16_t)half;
#endif
}

template <typename T> struct RowStore;

template <> struct RowStore<float> {
  static inline float put(float* row, int i, float x) { row[i] = x; return x; }
  static inline float get(const float* row, int i) { return row[i]; }
  static inline float dot(const float* v, const float* s, int k) {
    // 8 float lanes x 2 accumulators, summed in double (the order of the original AVX2 kernel)
    float a[16] = {0};
    int i = 0;
    for (; i + 16 <= k; i += 16)
      for (int l = 0; l < 16; ++l) a[l] += v[i + l] * s[i + l];
    for (; i + 8 <= k; i += 8)
      for (int l = 0; l < 8; ++l) a[l] += v[i + l] * s[i + l];
    double acc = 0;
    for (int l = 0; l < 8; ++l) acc += (double)(a[l] + a[l + 8]);
    for (; i < k; ++i) acc += (double)v[i] * s[i];
    return (float)acc;
  }
  static inline void axpy(float* s, const float* v, float a, int k) {
    for (int i = 0; i < k; ++i) s[i] += a * v[i];
  }
};

template <> struct RowStore<uint16_t> {
  static inline float put(uint16_t* row, int i, float x) {
    const uint16_t h = float_to_half(x);
    row[i] = h;
    return half_to_float(h);
  }
  static inline float get(const uint16_t* row, int i) { return half_to_float(row[i]); }
  static inline float dot(const uint16_t* v, const float* s, int k) {
    float a[16] = {0};
    int i = 0;
    for (; i + 16 <= k; i += 16)
      for (int l = 0; l < 16; ++l) a[l] += half_to_float(v[i + l]) * s[i + l];
    for (; i + 8 <= k; i += 8)
      for (int l = 0; l < 8; ++l) a[l] += half_to_float(v[i + l]) * s[i + l];
    double acc = 0;
    for (int l = 0; l < 8; ++l) acc += (double)(a[l] + a[l + 8]);
    for (; i < k; ++i) acc += (double)half_to_float(v[i]) * s[i];
    return (float)acc;
  }
  static inline void axpy(float* s, const uint16_t* v, float a, int k) {
    for (int i = 0; i < k; ++i) s[i] += a * half_to_float(v[i]);
  }
};

}  // namespace plink4r
#endif
