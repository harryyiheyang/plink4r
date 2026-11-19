// Products of a geno object (X = H + D, see geno_store.h): the kernels write into caller buffers
// (column-major); the R entry points below and the C-callable ones (geno_api.cpp) both call them.
#include "geno_store.h"
#include "geno_products.h"
#include "geno_kernel.h"

using namespace plink4r::gstore;

namespace {

constexpr std::size_t kDRows = 2048;  // rows per block of the D accumulation (a multiple of 64)
constexpr std::size_t kDCols = 256;   // D columns per work item of the D accumulation
constexpr std::size_t kDLine = 16;    // columns decoded together (one cache line of a row of Y)
constexpr std::size_t kDSpan = 1024;  // output rows per work item of the D accumulation

// s[a] += sum_k d[k] Y[r[k] * ld + a] for a < len, the terms added in the order of k (4 rows per pass over s).
using DAxpy = void (*)(const float*, std::size_t, const std::uint32_t*, const double*, std::size_t, double*,
                       std::size_t);

__attribute__((always_inline)) inline void d_axpy_body(const float* Y, std::size_t ld, const std::uint32_t* r,
                                                       const double* d, std::size_t m, double* s, std::size_t len) {
  std::size_t k = 0;
  for (; k + 4 <= m; k += 4) {
    const float *y0 = Y + r[k] * ld, *y1 = Y + r[k + 1] * ld, *y2 = Y + r[k + 2] * ld, *y3 = Y + r[k + 3] * ld;
    const double d0 = d[k], d1 = d[k + 1], d2 = d[k + 2], d3 = d[k + 3];
    #pragma omp simd
    for (std::size_t a = 0; a < len; ++a) {
      double v = s[a];
      v += d0 * y0[a];
      v += d1 * y1[a];
      v += d2 * y2[a];
      v += d3 * y3[a];
      s[a] = v;
    }
  }
  for (; k < m; ++k) {
    const float* y = Y + r[k] * ld;
    const double dk = d[k];
    #pragma omp simd
    for (std::size_t a = 0; a < len; ++a) s[a] += dk * y[a];
  }
}

void d_axpy_generic(const float* Y, std::size_t ld, const std::uint32_t* r, const double* d, std::size_t m,
                    double* s, std::size_t len) {
  d_axpy_body(Y, ld, r, d, m, s, len);
}

#ifdef GENO_KERNEL_X86
__attribute__((target("avx2,fma")))
void d_axpy_avx2(const float* Y, std::size_t ld, const std::uint32_t* r, const double* d, std::size_t m,
                 double* s, std::size_t len) {
  d_axpy_body(Y, ld, r, d, m, s, len);
}
#endif

// Chosen as select_kernel(): the AVX2 + FMA loop whenever la() picked avx2 or avx512 (an AVX-512 loop
// measured no faster: the loop is bound by reading Y).
DAxpy select_d_axpy() {
#ifdef GENO_KERNEL_X86
  if (std::strcmp(plink4r::la().isa, "generic") != 0) return d_axpy_avx2;
#endif
  return d_axpy_generic;
}

// X'X += H'D + D'H + D'D, as S + S' with S[a, b] = sum_i D_bi (H_ai + D_ai / 2). With cols (0-based, length
// pc), the columns cols of g in that order (out is pc x pc); without, all columns of g. With per-row weights
// wrow, D_bi is multiplied by wrow_i, giving X'WX += H'WD + D'WH + D'WD. Per block of kDRows rows, Y = (H + D / 2)
// / sc of every column (sc = unit / 2) is decoded once into a row-major float buffer (row i contiguous over the
// columns); Y = 2 H / unit + q is an integer below 2^17 for a genotype column and q for an appended one (H = 0),
// so float holds it exactly. Then work items of kDSpan rows x kDCols columns of S, shared out by OpenMP
// (disjoint outputs), add D_bi Y[i, ] to S[, b] for the nonzeros of their columns in the block, in double with
// the span of S and the rows of Y it reads in cache; sc is applied at the end. Every S[a, b] sums its terms in
// ascending row order.
void add_d_xtx(const Geno& g, double* out, int threads, const int* cols = nullptr, std::size_t pc = 0,
               const double* wrow = nullptr) {
  const std::size_t n = g.n, p = cols ? pc : g.p;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  if (!g.dnnz()) return;
  const DAxpy axpy = select_d_axpy();
  const std::size_t RB = std::min(kDRows, n);
  std::vector<double> S(p * p, 0.0), sc(p), val;  // val, row: the nonzeros of the block, column by column
  std::vector<float> Y(RB * p);
  std::vector<std::uint32_t> row;
  for (std::size_t a = 0; a < p; ++a) sc[a] = 0.5 * g.D[col(a)].unit;
  std::vector<std::size_t> e0(p), off(p + 1, 0);  // per column: first nonzero of the block, start in row / val
  const std::size_t ns = (p + kDSpan - 1) / kDSpan, ng = (p + kDCols - 1) / kDCols, nc = (p + kDLine - 1) / kDLine;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> yb(kDLine * RB);
    for (std::size_t i0 = 0; i0 < n; i0 += RB) {
      const std::size_t R = std::min(RB, n - i0), w0 = i0 / 64, w1 = (i0 + R + 63) / 64;
      #pragma omp for schedule(static)
      for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(p); ++aa) {
        const std::size_t a = static_cast<std::size_t>(aa);
        const DCol& c = g.D[col(a)];
        e0[a] = c.mask.empty() ? 0 : d_rank(c, w0);
        off[a + 1] = c.mask.empty() ? 0 : d_rank(c, w1) - e0[a];
      }
      #pragma omp single
      {
        std::partial_sum(off.begin(), off.end(), off.begin());
        if (row.size() < off[p]) { row.resize(off[p]); val.resize(off[p]); }
      }
      const bool any = off[p] > 0;
      #pragma omp barrier
      if (!any) continue;  // every thread has read off[p] before the next block rewrites it
      #pragma omp for schedule(static)
      for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(nc); ++cc) {
        const std::size_t a0 = static_cast<std::size_t>(cc) * kDLine, na = std::min(kDLine, p - a0);
        for (std::size_t t = 0; t < na; ++t) {
          const std::size_t a = a0 + t, j = col(a);
          double* y = yb.data() + t * RB;
          std::fill(y, y + R, 0.0);
          const u64* lo = g.G.lo(j);
          const u64* hi = g.G.hi(j);
          for (std::size_t w = w0; w < w1; ++w) {
            for (u64 bits = lo[w]; bits; bits &= bits - 1) y[64 * (w - w0) + __builtin_ctzll(bits)] = 1.0;
            for (u64 bits = hi[w]; bits; bits &= bits - 1) y[64 * (w - w0) + __builtin_ctzll(bits)] = 2.0;
          }
          const DCol& c = g.D[j];
          const double inv = 1.0 / sc[a];
          for (std::size_t r = 0; r < R; ++r) y[r] *= inv;
          std::size_t e = e0[a], k = off[a];
          for (std::size_t w = w0; w < w1 && off[a + 1] > off[a]; ++w)
            for (u64 bits = c.mask[w]; bits; bits &= bits - 1) {
              const std::size_t r = 64 * (w - w0) + __builtin_ctzll(bits);
              const int q = c.q[e++];
              row[k] = static_cast<std::uint32_t>(r);
              val[k++] = wrow ? q * c.unit * wrow[i0 + r] : q * c.unit;
              y[r] += q;
            }
        }
        for (std::size_t r = 0; r < R; ++r)
          for (std::size_t t = 0; t < na; ++t) Y[r * p + a0 + t] = static_cast<float>(yb[t * RB + r]);
      }
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t it = 0; it < static_cast<std::ptrdiff_t>(ns * ng); ++it) {
        const std::size_t s0 = (static_cast<std::size_t>(it) % ns) * kDSpan, len = std::min(kDSpan, p - s0);
        const std::size_t b0 = (static_cast<std::size_t>(it) / ns) * kDCols, b1 = std::min(p, b0 + kDCols);
        for (std::size_t b = b0; b < b1; ++b)
          if (off[b + 1] > off[b])
            axpy(Y.data() + s0, p, row.data() + off[b], val.data() + off[b], off[b + 1] - off[b], S.data() + s0 + b * p,
                 len);
      }  // implicit barrier: Y and the lists are rewritten for the next block after every work item of this one
    }
    #pragma omp for schedule(static)
    for (std::ptrdiff_t bb = 0; bb < static_cast<std::ptrdiff_t>(p); ++bb) {
      const std::size_t b = static_cast<std::size_t>(bb);
      for (std::size_t a = 0; a < p; ++a) out[a + b * p] += sc[a] * S[a + b * p] + sc[b] * S[b + a * p];
    }
  }
}

}  // namespace

namespace {

// The kernels of the products, writing into caller buffers (column-major); the R functions below and the
// C-callable entry points (bottom of this file) both call them.

// Rows r0..r0+R-1 of stored column j (H + D) as T (float or double) into dst[0..R). grid: the value in steps of
// 0.01 minus 100 (an integer in -100..100), for a column on the 0.01 grid (on_grid).
template <class T>
void decode_rows(const Geno& g, std::size_t j, std::size_t r0, std::size_t R, T* dst, bool grid = false) {
  const T v0 = grid ? T(-kDosSteps) : T(0), v1 = grid ? T(0) : T(1), v2 = grid ? T(kDosSteps) : T(2);
  std::fill(dst, dst + R, v0);
  const u64* lo = g.G.lo(j);
  const u64* hi = g.G.hi(j);
  const std::size_t w0 = r0 / 64, w1 = (r0 + R + 63) / 64;
  for (std::size_t x = w0; x < w1; ++x) {
    for (u64 bits = lo[x]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * x + __builtin_ctzll(bits);
      if (i >= r0 && i < r0 + R) dst[i - r0] = v1;
    }
    for (u64 bits = hi[x]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * x + __builtin_ctzll(bits);
      if (i >= r0 && i < r0 + R) dst[i - r0] = v2;
    }
  }
  const DCol& dc = g.D[j];
  if (!dc.mask.empty()) {
    std::size_t e = d_rank(dc, w0);
    for (std::size_t x = w0; x < w1 && x < dc.mask.size(); ++x)
      for (u64 bits = dc.mask[x]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * x + __builtin_ctzll(bits);
        const double d = grid ? dc.q[e++] : dc.q[e++] * dc.unit;
        if (i >= r0 && i < r0 + R) dst[i - r0] += static_cast<T>(d);
      }
  }
}

// The dense product loop of xtwx_float_into and xtx_dense_into, into o (P x P; P = pc for columns cols, else
// g.p). Per block of RB rows, every column is decoded once into a shared buffer of T (column a at X[a * LD],
// LD = RB + pad: no power-of-two stride, which thrashes the cache); then the output tiles of TS x TS of the
// lower triangle are shared out by OpenMP (disjoint outputs), each adding its block product to o through
// tile(X, LD, r0, R, a0, na, b0, nb) of the thread's own worker from make(). The lower triangle is mirrored.
// grid: the columns are decoded as centred 0.01 steps (decode_rows).
template <class T, class Make>
void dense_tiles(const Geno& g, const int* cols, std::size_t pc, std::size_t RB, std::size_t pad, std::size_t TS,
                 double* o, int threads, Make make, bool grid = false) {
  const std::size_t n = g.n, P = cols ? pc : g.p, LD = RB + pad;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  const std::size_t nt = (P + TS - 1) / TS;
  std::vector<std::pair<std::size_t, std::size_t>> tiles;
  for (std::size_t b = 0; b < nt; ++b) for (std::size_t a = b; a < nt; ++a) tiles.emplace_back(a, b);
  std::vector<T> X(LD * P);
  std::fill(o, o + P * P, 0.0);
  #pragma omp parallel num_threads(threads)
  {
    auto tile = make();
    for (std::size_t r0 = 0; r0 < n; r0 += RB) {
      const std::size_t R = std::min(RB, n - r0);
      #pragma omp for schedule(static)
      for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(P); ++aa)
        decode_rows(g, col(static_cast<std::size_t>(aa)), r0, R, &X[static_cast<std::size_t>(aa) * LD], grid);
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t tt = 0; tt < static_cast<std::ptrdiff_t>(tiles.size()); ++tt) {
        const std::size_t a0 = tiles[tt].first * TS, b0 = tiles[tt].second * TS;
        tile(X.data(), LD, r0, R, a0, std::min(TS, P - a0), b0, std::min(TS, P - b0));
      }  // implicit barrier: the buffer is rewritten for the next block after every tile of this one
    }
  }
  for (std::size_t b = 0; b < P; ++b)  // mirror the lower triangle
    for (std::size_t a = b + 1; a < P; ++a) o[b + a * P] = o[a + b * P];
}

constexpr std::size_t kFTile = 256;    // columns per side of an output tile of the float kernel
constexpr std::size_t kFRows = 4096;   // rows per block of the float kernel
constexpr std::size_t kFPad = 16;      // floats after each column of the shared buffer (no power-of-two stride)

}  // namespace

namespace plink4r {
namespace gstore {

// X'WX into o (p x p, or pc x pc for columns cols), W = diag(w), any real w, in float products with double
// sums (dense_tiles). Blocks of kFRows rows, 16 KB per column of the shared buffer (66 MB at p = 4000, 329 MB
// at p = 20000); a tile copies its row group transposed and times w, multiplies it with its column group in
// the shared buffer by the run-time-dispatched serial sgemm and accumulates in double.
void xtwx_float_into(const Geno& g, const double* w, const int* cols, std::size_t pc, double* o, int threads) {
  const std::size_t P = cols ? pc : g.p, RB = std::min(kFRows, g.n);
  const plink4r::LaFns& LA = plink4r::la();
  dense_tiles<float>(g, cols, pc, RB, kFPad, kFTile, o, threads, [&] {
    return [&, At = std::vector<float>(kFTile * RB), Pt = std::vector<float>(kFTile * kFTile)](
               const float* X, std::size_t LD, std::size_t r0, std::size_t R, std::size_t a0, std::size_t na,
               std::size_t b0, std::size_t nb) mutable {
      const float* xa = &X[a0 * LD];
      for (std::size_t r = 0; r < R; ++r) {  // A' : na x R, ld kFTile
        const double wr = w[r0 + r];
        for (std::size_t c = 0; c < na; ++c) At[c + r * kFTile] = static_cast<float>(xa[c * LD + r] * wr);
      }
      LA.sgemm_nn(At.data(), static_cast<int>(na), static_cast<int>(R), static_cast<int>(kFTile), &X[b0 * LD],
                  static_cast<int>(nb), static_cast<int>(LD), Pt.data(), static_cast<int>(kFTile));
      for (std::size_t y = 0; y < nb; ++y)
        for (std::size_t x = 0; x < na; ++x) o[(a0 + x) + (b0 + y) * P] += Pt[x + y * kFTile];
    };
  });
}

}  // namespace gstore
}  // namespace plink4r

namespace {


constexpr std::size_t kXTile = 256;    // columns per side of an output tile of the double X'X kernel
constexpr std::size_t kXRows = 2048;   // rows per block of the double X'X kernel
constexpr std::size_t kXPad = 8;       // doubles after each column of the shared buffer (no power-of-two stride)

// X'X into o (p x p, or pc x pc for columns cols) in double (dense_tiles): blocks of kXRows rows, 16 KB per
// column of the shared buffer (66 MB at p = 4000, 329 MB at p = 20000); a tile adds X_a'X_b of its row and
// column groups to o by the run-time-dispatched serial dgemm. Columns in 1/16384 steps hold multiples of 2^-14
// in [0, 2], so for n below 2^23 every product and partial sum is exact and the result equals the popcount one.
void xtx_dense_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads) {
  const std::size_t P = cols ? pc : g.p, RB = std::min(kXRows, g.n);
  const plink4r::LaFns& LA = plink4r::la();
  dense_tiles<double>(g, cols, pc, RB, kXPad, kXTile, o, threads, [&] {
    return [&](const double* X, std::size_t LD, std::size_t, std::size_t R, std::size_t a0, std::size_t na,
               std::size_t b0, std::size_t nb) {
      LA.dgemm_tn(&X[a0 * LD], static_cast<int>(R), static_cast<int>(na), static_cast<int>(LD), &X[b0 * LD],
                  static_cast<int>(nb), static_cast<int>(LD), &o[a0 + b0 * P], static_cast<int>(P), 1);
    };
  });
}

constexpr double kDDenseFrac = 0.01;     // D nonzero fraction above which xtwx_bin_into uses the float kernel
constexpr std::size_t kGRows = 1024;   // rows per block of the grid X'X kernel: partial sums <= 1024 * 100^2 < 2^24

// Whether every column of cols (or all) is on the 0.01 grid: no D, or D in steps of 0.01 (dosage PGEN).
bool on_grid(const Geno& g, const int* cols, std::size_t pc) {
  const std::size_t P = cols ? pc : g.p;
  for (std::size_t a = 0; a < P; ++a) {
    const DCol& c = g.D[cols ? cols[a] : a];
    if (!c.mask.empty() && c.unit != kDosUnit) return false;
  }
  return true;
}

// X'X into o (p x p, or pc x pc for columns cols) for columns on the 0.01 grid (on_grid), exact, in float
// products (dense_tiles, tiles and buffer as xtwx_float_into, blocks of kGRows rows). A column is decoded as
// the integer c = 100 x - 100 in -100..100; within a block every product and partial sum is an integer below
// 2^24, so the serial sgemm is exact in any order; the blocks are summed in double (integers below 2^53) and
// X'X = (C'C + 100 (s_a + s_b) - 1e4 n) / 1e4, s the column sums of 100 x.
void xtx_grid_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads) {
  const std::size_t P = cols ? pc : g.p, RB = std::min(kGRows, g.n);
  const plink4r::LaFns& LA = plink4r::la();
  dense_tiles<float>(g, cols, pc, RB, kFPad, kFTile, o, threads, [&] {
    return [&, At = std::vector<float>(kFTile * RB), Pt = std::vector<float>(kFTile * kFTile)](
               const float* X, std::size_t LD, std::size_t, std::size_t R, std::size_t a0, std::size_t na,
               std::size_t b0, std::size_t nb) mutable {
      const float* xa = &X[a0 * LD];
      for (std::size_t r = 0; r < R; ++r)  // A' : na x R, ld kFTile
        for (std::size_t c = 0; c < na; ++c) At[c + r * kFTile] = xa[c * LD + r];
      LA.sgemm_nn(At.data(), static_cast<int>(na), static_cast<int>(R), static_cast<int>(kFTile), &X[b0 * LD],
                  static_cast<int>(nb), static_cast<int>(LD), Pt.data(), static_cast<int>(kFTile));
      for (std::size_t y = 0; y < nb; ++y)
        for (std::size_t x = 0; x < na; ++x) o[(a0 + x) + (b0 + y) * P] += Pt[x + y * kFTile];
    };
  }, true);
  std::vector<double> s(P);
  for (std::size_t a = 0; a < P; ++a)
    s[a] = static_cast<double>(step_sum_col(g, cols ? static_cast<std::size_t>(cols[a]) : a, kDosSteps));
  const double k = kDosSteps, nk2 = k * k * static_cast<double>(g.n);
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t b = 0; b < static_cast<std::ptrdiff_t>(P); ++b)
    for (std::size_t a = 0; a < P; ++a) {
      double& v = o[a + static_cast<std::size_t>(b) * P];
      v = (v + k * (s[a] + s[static_cast<std::size_t>(b)]) - nk2) / (k * k);
    }
}

constexpr double kDDenseFracXtx = 0.07;  // D nonzero fraction above which xtx_into uses a dense kernel

// Whether more than frac of the cells of the columns cols (or all) have a nonzero D.
bool d_dense(const Geno& g, const int* cols, std::size_t pc, double frac) {
  const std::size_t P = cols ? pc : g.p;
  double nnz = 0.0;
  for (std::size_t a = 0; a < P; ++a) nnz += static_cast<double>(g.D[cols ? cols[a] : a].q.size());
  return nnz > frac * static_cast<double>(g.n) * static_cast<double>(P);
}

}  // namespace

namespace plink4r {
namespace gstore {

// X'X of the stored matrix into o (p x p), or of columns cols (0-based, length pc) into o (pc x pc): exact
// popcounts of H'H plus add_d_xtx, or, when more than kDDenseFracXtx of the cells of the columns have a nonzero
// D (dosages), a dense kernel: xtx_grid_into for columns on the 0.01 grid, else the double xtx_dense_into.
void xtx_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads) {
  if (d_dense(g, cols, pc, kDDenseFracXtx)) {
    if (on_grid(g, cols, pc)) xtx_grid_into(g, cols, pc, o, threads);
    else xtx_dense_into(g, cols, pc, o, threads);
    return;
  }
  const TileKernel kernel = select_kernel();
  if (!cols) {
    const std::size_t p = g.p;
    dot_block(g.G, 0, p, g.G, 0, p, true, o, p, threads, kernel);
    add_d_xtx(g, o, threads);
    return;
  }
  Planes sub;
  sub.words = g.G.words;
  sub.bits.resize(2 * pc * sub.words);
  for (std::size_t a = 0; a < pc; ++a) {
    const std::size_t j = static_cast<std::size_t>(cols[a]);
    std::copy(g.G.lo(j), g.G.lo(j) + sub.words, sub.bits.begin() + 2 * a * sub.words);
    std::copy(g.G.hi(j), g.G.hi(j) + sub.words, sub.bits.begin() + (2 * a + 1) * sub.words);
  }
  dot_block(sub, 0, pc, sub, 0, pc, true, o, pc, threads, kernel);
  add_d_xtx(g, o, threads, cols, pc);
}

}  // namespace gstore
}  // namespace plink4r

namespace {


constexpr std::size_t kBinRows = 64;     // rows per weight bin: one 64-bit word of the gathered planes
constexpr std::size_t kBinGroup = 384;   // bins per streamed group (24576 rows; a multiple of 4)
}  // namespace

namespace plink4r {
namespace gstore {

// X'W~X into o (p x p, or pc x pc for columns cols), W~ = diag(w binned): rows sorted by w (ties by row
// index), each run of kBinRows sorted rows takes the mean w of its rows. With the weight constant within a
// bin, H'W~H = sum_b wbar_b H_b'H_b by exact popcounts: per group of kBinGroup bins the rows of the group
// are gathered for every column (one bin per word), then the 64 x 64 tiles of the upper triangle are shared
// out by OpenMP (disjoint outputs). The D terms come from add_d_xtx with the binned per-row weights. When
// more than kDDenseFrac of the cells of the columns have a nonzero D (dosages), the popcounts do not pay:
// X'WX with w itself comes from xtwx_float_into, and W~ = W. w_bin (n, row order, may be null) receives the
// weights of the result.
void xtwx_bin_into(const Geno& g, const double* w, const int* cols, std::size_t pc, double* o, double* w_bin,
                   int threads) {
  const std::size_t n = g.n, P = cols ? pc : g.p;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  if (d_dense(g, cols, pc, kDDenseFrac)) {
    xtwx_float_into(g, w, cols, pc, o, threads);
    if (w_bin) std::copy(w, w + n, w_bin);
    return;
  }
  std::vector<std::uint32_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0u);
  std::stable_sort(ord.begin(), ord.end(), [w](std::uint32_t a, std::uint32_t b) { return w[a] < w[b]; });
  const std::size_t nbin = (n + kBinRows - 1) / kBinRows, nbin4 = (nbin + 3) / 4 * 4;
  std::vector<double> wbar(nbin4, 0.0), wq(n);  // padding bins (no rows) weigh 0
  for (std::size_t b = 0; b < nbin; ++b) {
    const std::size_t r0 = b * kBinRows, r1 = std::min(n, r0 + kBinRows);
    double s = 0.0;
    for (std::size_t r = r0; r < r1; ++r) s += w[ord[r]];
    wbar[b] = s / static_cast<double>(r1 - r0);
    for (std::size_t r = r0; r < r1; ++r) wq[ord[r]] = wbar[b];
  }
  if (w_bin) std::copy(wq.begin(), wq.end(), w_bin);
  std::fill(o, o + P * P, 0.0);
  const WTileKernel kernel = select_wkernel();
  const std::size_t nt = (P + kTile - 1) / kTile, gw = std::min(kBinGroup, nbin4);
  std::vector<u64> store(2 * P * gw + 3);  // 3 spare words to start the buffer on 32 bytes
  void* start = store.data();
  std::size_t space = store.size() * sizeof(u64);
  u64* buf = static_cast<u64*>(std::align(32, 2 * P * gw * sizeof(u64), start, space));
  for (std::size_t b0 = 0; b0 < nbin4; b0 += kBinGroup) {
    const std::size_t W = std::min(kBinGroup, nbin4 - b0);  // words per plane of this group
    const std::size_t g0 = b0 * kBinRows, g1 = std::min(n, (b0 + W) * kBinRows);
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(P); ++aa) {
      const std::size_t a = static_cast<std::size_t>(aa), j = col(a);
      const u64* lo = g.G.lo(j);
      const u64* hi = g.G.hi(j);
      u64* ol = buf + 2 * a * W;
      u64* oh = ol + W;
      for (std::size_t k = 0; k < W; ++k) {
        u64 l = 0, h = 0;
        const std::size_t r0 = g0 + kBinRows * k, r1 = std::min(g1, r0 + kBinRows);
        for (std::size_t r = r0; r < r1; ++r) {
          const std::uint32_t i = ord[r];
          l |= ((lo[i >> 6] >> (i & 63u)) & 1u) << (r - r0);
          h |= ((hi[i >> 6] >> (i & 63u)) & 1u) << (r - r0);
        }
        ol[k] = l;
        oh[k] = h;
      }
    }
    #pragma omp parallel num_threads(threads)
    {
      std::vector<double> acc(kTile * kTile);
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t t = 0; t < static_cast<std::ptrdiff_t>(nt * nt); ++t) {
        const std::size_t ti = static_cast<std::size_t>(t) / nt, tj = static_cast<std::size_t>(t) % nt;
        if (ti > tj) continue;
        const std::size_t i0 = ti * kTile, j0 = tj * kTile, ni = std::min(kTile, P - i0), nj = std::min(kTile, P - j0);
        std::fill(acc.begin(), acc.end(), 0.0);
        kernel(buf, W, i0, ni, j0, nj, wbar.data() + b0, acc.data());
        for (std::size_t x = 0; x < ni; ++x)
          for (std::size_t y = 0; y < nj; ++y) o[(i0 + x) + (j0 + y) * P] += acc[x * kTile + y];
      }
    }
  }
  for (std::size_t b = 0; b < P; ++b)  // mirror the upper triangle
    for (std::size_t a = b + 1; a < P; ++a) o[a + b * P] = o[b + a * P];
  add_d_xtx(g, o, threads, cols, pc, wq.data());
}

}  // namespace gstore
}  // namespace plink4r

namespace {


}  // namespace

namespace plink4r {
namespace gstore {

// Standardizes a weighted cross product in place: o (pc x pc, columns cols or all when null) becomes
// (X'WX - c s' - s c' + (sum w) c c') / (d d'), s = X'w in double from the planes and D; center c and div d
// are per output column (div may be null: no division).
void xtwx_std(const Geno& g, const double* w, const int* cols, std::size_t pc, const double* center,
              const double* div, double* o, int threads) {
  const std::size_t n = g.n, P = cols ? pc : g.p;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  std::vector<double> sw(P, 0.0);
  double wsum = 0.0;
  for (std::size_t i = 0; i < n; ++i) wsum += w[i];
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(P); ++aa) {
    const std::size_t a = static_cast<std::size_t>(aa), j = col(a);
    double acc = 0.0;
    const u64* lo = g.G.lo(j);
    const u64* hi = g.G.hi(j);
    for (std::size_t x = 0; x < g.G.words; ++x) {
      for (u64 bits = lo[x]; bits; bits &= bits - 1) acc += w[64 * x + __builtin_ctzll(bits)];
      for (u64 bits = hi[x]; bits; bits &= bits - 1) acc += 2.0 * w[64 * x + __builtin_ctzll(bits)];
    }
    const DCol& dc = g.D[j];
    std::size_t e = 0;
    for (std::size_t x = 0; x < dc.mask.size(); ++x)
      for (u64 bits = dc.mask[x]; bits; bits &= bits - 1) acc += dc.q[e++] * dc.unit * w[64 * x + __builtin_ctzll(bits)];
    sw[a] = acc;
  }
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t bb = 0; bb < static_cast<std::ptrdiff_t>(P); ++bb) {
    const std::size_t b = static_cast<std::size_t>(bb);
    for (std::size_t a = 0; a < P; ++a) {
      double v = o[a + b * P] - center[a] * sw[b] - sw[a] * center[b] + wsum * center[a] * center[b];
      if (div) v /= div[a] * div[b];
      o[a + b * P] = v;
    }
  }
}

}  // namespace gstore
}  // namespace plink4r

namespace {


}  // namespace

namespace plink4r {
namespace gstore {

// X'M into o (p x k) for a dense n x k matrix m.
void xtm_into(const Geno& g_, const double* m, std::size_t k, double* o, int threads) {
  const Geno* g = &g_;
  const std::size_t n = g->n, p = g->p, W = g->G.words;
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    std::vector<double> s1(k, 0.0), s2(k, 0.0), sd(k, 0.0);
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    for (std::size_t w = 0; w < W; ++w) {
      for (u64 bits = lo[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) s1[c] += m[i + c * n];
      }
      for (u64 bits = hi[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) s2[c] += m[i + c * n];
      }
    }
    const DCol& dc = g->D[j];
    std::size_t e = 0;
    for (std::size_t w = 0; w < dc.mask.size(); ++w)
      for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) {
        const double d = dc.q[e++] * dc.unit;
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) sd[c] += d * m[i + c * n];
      }
    for (std::size_t c = 0; c < k; ++c) o[j + c * p] = s1[c] + 2.0 * s2[c] + sd[c];
  }
}

}  // namespace gstore
}  // namespace plink4r

namespace {


}  // namespace

namespace plink4r {
namespace gstore {

// X B into o (n x k, overwritten) for a dense p x k matrix b; all-zero rows of b are skipped.
void mx_into(const Geno& g_, const double* b, std::size_t k, double* o, int threads) {
  const Geno* g = &g_;
  const std::size_t n = g->n, p = g->p;
  std::fill(o, o + n * k, 0.0);
  std::vector<char> nonzero(p, 0);
  for (std::size_t j = 0; j < p; ++j)
    for (std::size_t c = 0; c < k; ++c)
      if (b[j + c * p] != 0.0) { nonzero[j] = 1; break; }
  const std::size_t used_words = (n + 63) / 64;
  #pragma omp parallel num_threads(threads)
  {
    const int nt = plink4r::thread_team(), t = plink4r::thread_id();
    const std::size_t w0 = used_words * t / nt, w1 = used_words * (t + 1) / nt;
    std::vector<double> bj(k);
    for (std::size_t j = 0; j < p; ++j) {
      if (!nonzero[j]) continue;
      for (std::size_t c = 0; c < k; ++c) bj[c] = b[j + c * p];
      const u64* lo = g->G.lo(j);
      const u64* hi = g->G.hi(j);
      for (std::size_t w = w0; w < w1; ++w) {
        for (u64 bits = lo[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) o[i + c * n] += bj[c];
        }
        for (u64 bits = hi[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) o[i + c * n] += 2.0 * bj[c];
        }
      }
      const DCol& dc = g->D[j];
      if (!dc.mask.empty() && w0 < w1) {
        std::size_t e = d_rank(dc, w0);
        for (std::size_t w = w0; w < w1; ++w)
          for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) {
            const double d = dc.q[e++] * dc.unit;
            const std::size_t i = 64 * w + __builtin_ctzll(bits);
            for (std::size_t c = 0; c < k; ++c) o[i + c * n] += d * bj[c];
          }
      }
    }
  }
}

}  // namespace gstore
}  // namespace plink4r

namespace {


}  // namespace

namespace plink4r {
namespace gstore {

// Dense block X[rows, cols] of the stored matrix (0-based) into o (nr x nc).
void dense_into(const Geno& g_, const int* rows, std::size_t nr, const int* cols, std::size_t nc, double* o,
                int threads) {
  const Geno* g = &g_;
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(nc); ++cc) {
    const std::size_t j = static_cast<std::size_t>(cols[cc]);
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    double* col = o + static_cast<std::size_t>(cc) * nr;
    for (std::size_t r = 0; r < nr; ++r) {
      const std::size_t i = static_cast<std::size_t>(rows[r]);
      const std::size_t w = i / 64;
      const u64 bit = u64(1) << (i % 64);
      col[r] = (lo[w] & bit) ? 1.0 : (hi[w] & bit) ? 2.0 : 0.0;
    }
    const DCol& dc = g->D[j];
    if (!dc.mask.empty())
      for (std::size_t r = 0; r < nr; ++r) {
        const double d = d_at(dc, static_cast<std::size_t>(rows[r]));
        if (d != 0.0) col[r] += d;
      }
  }
}

}  // namespace gstore
}  // namespace plink4r

namespace {


}  // namespace

namespace plink4r {
namespace gstore {

// o = (k - s u v') / (da db'), elementwise (r x c); a null da or db is 1. o may be k (in place).
void std_into(const double* k, std::size_t r, std::size_t c, const double* u, const double* v, double s,
              const double* da, const double* db, double* o, int threads) {
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(c); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    const double vj = s * v[j], bj = db ? db[j] : 1.0;
    for (std::size_t i = 0; i < r; ++i) o[i + j * r] = (k[i + j * r] - u[i] * vj) / ((da ? da[i] : 1.0) * bj);
  }
}

}  // namespace gstore
}  // namespace plink4r

namespace {



}  // namespace

// X'X of the stored, unstandardized matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtx_cpp(SEXP ptr, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t p = g->p;
  Rcpp::NumericMatrix out(p, p);
  xtx_into(*g, nullptr, 0, out.begin(), threads);
  return out;
}

// X'M for a dense n x k matrix M.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtm_cpp(SEXP ptr, const Rcpp::NumericMatrix& M, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  const std::size_t k = M.ncol();
  if (static_cast<std::size_t>(M.nrow()) != n) Rcpp::stop("M must have one row per sample.");
  Rcpp::NumericMatrix out(p, k);
  xtm_into(*g, M.begin(), k, out.begin(), threads);
  return out;
}

// X B for a dense p x k matrix B; all-zero rows of B are skipped.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_mx_cpp(SEXP ptr, const Rcpp::NumericMatrix& B, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  const std::size_t k = B.ncol();
  if (static_cast<std::size_t>(B.nrow()) != p) Rcpp::stop("B must have one row per variant.");
  Rcpp::NumericMatrix out(n, k);
  mx_into(*g, B.begin(), k, out.begin(), threads);
  return out;
}

// Raw A'B of two geno objects over the same samples, both without D (hard
// calls, median fill): exact popcount on the planes.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_cross_cpp(SEXP ptr_a, SEXP ptr_b, int threads) {
  Rcpp::XPtr<Geno> a(ptr_a), b(ptr_b);
  threads = resolve_threads(threads);
  if (a->n != b->n) Rcpp::stop("A and B must have the same samples.");
  if (a->dnnz() || b->dnnz()) Rcpp::stop("geno_cross_cpp needs hard-call genotypes without D.");
  const TileKernel kernel = select_kernel();
  Rcpp::NumericMatrix out(a->p, b->p);
  dot_block(a->G, 0, a->p, b->G, 0, b->p, false, out.begin(), a->p, threads, kernel);
  return out;
}

// Dense block X[rows, cols] of the stored, unstandardized matrix (0-based).
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_dense_cpp(SEXP ptr, const Rcpp::IntegerVector& rows,
                                   const Rcpp::IntegerVector& cols, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t nr = rows.size(), nc = cols.size();
  Rcpp::NumericMatrix out(nr, nc);
  dense_into(*g, rows.begin(), nr, cols.begin(), nc, out.begin(), threads);
  return out;
}

// Standardized products in C++ (no BLAS from R): out = (K - s u v') / (da db'), elementwise; an empty
// da or db is 1. geno_crossprod / geno_cor give their centered, scaled forms through this.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_std_cpp(const Rcpp::NumericMatrix& K, const Rcpp::NumericVector& u,
                                 const Rcpp::NumericVector& v, double s, const Rcpp::NumericVector& da,
                                 const Rcpp::NumericVector& db, int threads) {
  threads = resolve_threads(threads);
  const std::size_t r = K.nrow(), c = K.ncol();
  if (static_cast<std::size_t>(u.size()) != r || static_cast<std::size_t>(v.size()) != c) Rcpp::stop("u, v do not match K.");
  const bool ha = da.size() > 0, hb = db.size() > 0;
  if ((ha && static_cast<std::size_t>(da.size()) != r) || (hb && static_cast<std::size_t>(db.size()) != c)) Rcpp::stop("da, db do not match K.");
  Rcpp::NumericMatrix out(r, c);
  std_into(K.begin(), r, c, u.begin(), v.begin(), s, ha ? da.begin() : nullptr, hb ? db.begin() : nullptr, out.begin(),
           threads);
  return out;
}

// A'A of a dense n x k matrix (double, columns pairs shared out by OpenMP).
// [[Rcpp::export]]
Rcpp::NumericMatrix dense_xtx_cpp(const Rcpp::NumericMatrix& A, int threads) {
  threads = resolve_threads(threads);
  const std::size_t n = A.nrow(), k = A.ncol();
  Rcpp::NumericMatrix out(k, k);
  const double* a = A.begin();
  double* o = out.begin();
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(k); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    for (std::size_t i = j; i < k; ++i) {
      double acc = 0.0;
      for (std::size_t r = 0; r < n; ++r) acc += a[r + i * n] * a[r + j * n];
      o[i + j * k] = acc; o[j + i * k] = acc;
    }
  }
  return out;
}

// X'WX of the stored matrix, W = diag(w), any real w (xtwx_float_into). With center (and div) given, returns
// the product of the standardized matrix: (X'WX - c s' - s c' + (sum w) c c') / (d d'), s = X'w.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtwx_cpp(SEXP ptr, const Rcpp::NumericVector& w, const Rcpp::NumericVector& center,
                                  const Rcpp::NumericVector& div, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  if (static_cast<std::size_t>(w.size()) != n) Rcpp::stop("w must have one value per sample.");
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(w[i])) Rcpp::stop("w must be finite.");
  if (center.size() > 0 && static_cast<std::size_t>(center.size()) != p) Rcpp::stop("center must have one value per column.");
  Rcpp::NumericMatrix out(p, p);
  xtwx_float_into(*g, w.begin(), nullptr, p, out.begin(), threads);
  if (center.size() > 0)
    xtwx_std(*g, w.begin(), nullptr, p, center.begin(), div.size() > 0 ? div.begin() : nullptr, out.begin(), threads);
  return out;
}

// X'W~X of the stored matrix with w binned in groups of 64 rows, or with w itself when more than 1% of the
// cells have a nonzero D (xtwx_bin_into), standardized as in geno_xtwx_cpp with those weights when center is
// given; w = the weights of the result in row order.
// [[Rcpp::export]]
Rcpp::List geno_xtwx_bin_cpp(SEXP ptr, const Rcpp::NumericVector& w, const Rcpp::NumericVector& center,
                             const Rcpp::NumericVector& div, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  if (static_cast<std::size_t>(w.size()) != n) Rcpp::stop("w must have one value per sample.");
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(w[i])) Rcpp::stop("w must be finite.");
  if (center.size() > 0 && static_cast<std::size_t>(center.size()) != p) Rcpp::stop("center must have one value per column.");
  Rcpp::NumericMatrix out(p, p);
  Rcpp::NumericVector wq(n);
  xtwx_bin_into(*g, w.begin(), nullptr, p, out.begin(), wq.begin(), threads);
  if (center.size() > 0)
    xtwx_std(*g, wq.begin(), nullptr, p, center.begin(), div.size() > 0 ? div.begin() : nullptr, out.begin(), threads);
  return Rcpp::List::create(Rcpp::Named("XtWX") = out, Rcpp::Named("w") = wq);
}
