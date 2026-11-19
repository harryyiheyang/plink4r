// IEEE binary16 <-> float, shared by vi_store.h and the row kernels in la_kernel.cpp. Static (internal linkage) so
// that each ISA copy of la_kernel.cpp keeps its own: F16C when the translation unit targets it, otherwise the exact
// scalar conversion (round to nearest even); both give the same values.
#ifndef PLINK4R_HALF_H
#define PLINK4R_HALF_H

#include <cstdint>
#include <cstring>
#if defined(__F16C__)
#include <immintrin.h>
#endif

namespace plink4r {

static inline float half_to_float(uint16_t h) {
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

static inline uint16_t float_to_half(float x) {
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

}  // namespace plink4r
#endif
