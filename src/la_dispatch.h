// Run-time choice of the dense linear-algebra kernels in la_kernel.cpp (compiled for generic, AVX2+FMA and
// AVX-512 by src/Makevars). la() returns the best table the CPU supports; PLINK4R_SIMD=avx512|avx2|scalar
// caps the choice (the same variable as the LD plane-dot kernel).
#ifndef PLINK4R_LA_DISPATCH_H
#define PLINK4R_LA_DISPATCH_H

#include <cstdlib>
#include <string>

namespace plink4r {

struct LaFns {
  int (*build_C)(const float* U, int ld_m, int k, const float* lambda, double delta, const int* rows, int nrows,
                 float* C, int threads, char* err);
  void (*sgemm_nn)(const float* A, int M, int K, int lda, const float* B, int N, int ldb, float* Out, int ldo);
  int (*sumstat_block)(const float* U, int m, int k, const float* lambda, const int* typed, int nt, const float* z,
                       const int* miss, int nm, float delta, float* zout, float* r2out, char* err);
  const char* isa;
};

}  // namespace plink4r

extern const plink4r::LaFns plink4r_la_generic;
#if defined(PLINK4R_LA_X86)
extern const plink4r::LaFns plink4r_la_avx2;
extern const plink4r::LaFns plink4r_la_avx512;
#endif

namespace plink4r {

inline const LaFns& pick_la() {
  const char* env = std::getenv("PLINK4R_SIMD");
  const std::string want = env ? env : "";
#if defined(PLINK4R_LA_X86) && (defined(__GNUC__) || defined(__clang__))
  __builtin_cpu_init();
  const bool a512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
                    __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw");
  const bool a2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
  if (a512 && (want.empty() || want == "avx512")) return plink4r_la_avx512;
  if (a2 && (want.empty() || want == "avx512" || want == "avx2")) return plink4r_la_avx2;
#endif
  return plink4r_la_generic;
}
inline const LaFns& la() { static const LaFns& t = pick_la(); return t; }

}  // namespace plink4r
#endif
