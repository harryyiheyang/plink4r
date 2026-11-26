// The reading layer of the association scans (assoc_stream.cpp, gwas_kin_stream.cpp, gxe_kin_stream.cpp,
// xqtl_stream.cpp): the genotype file of a scan and its variants, the genotypes of one thread's block of variants in
// compact form, and the block size that fits the memory budget. Each scan computes its own statistics from the centred
// genotypes this layer hands out.
#ifndef PLINK4R_SCAN_BLOCK_H
#define PLINK4R_SCAN_BLOCK_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include "pgen_reader.h"
#include "variant_text.h"

namespace plink4r {

constexpr int kScanRowChunk = 2048;   // rows of G expanded at a time (multiple of 4)
constexpr int kScanBlockMax = 4096;   // largest automatic block per thread

// The genotype file of a scan: a BED (A1 = .bim column 5 counted) or a PGEN (ALT counted; hard calls or dosages), and
// the file rows (0-based) of the variants to scan and of the analysis samples.
struct ScanFile {
  std::unique_ptr<PgenFile> pf;
  std::unique_ptr<BedFile> bf;
  std::vector<int> vi, si;
  int n_raw = 0, ns = 0;
  bool ident = false, dosage = false;
  // m: the variants of the file. allele_ct: the number of alleles of every variant of a PGEN, by file row (a BED
  // has two of each; not read). var_rows: the file rows of the variants to scan. T: threads that will read.
  ScanFile(const std::string& path, bool bed, int n_raw_, std::size_t m, const std::vector<int>& allele_ct,
           std::vector<int> var_rows, const Rcpp::IntegerVector& sample_idx, int T)
      : vi(std::move(var_rows)), si(sample_idx.begin(), sample_idx.end()), n_raw(n_raw_),
        ns(static_cast<int>(sample_idx.size())) {
    ident = identity_samples(si.data(), ns, n_raw);
    if (bed) bf.reset(new BedFile(path, n_raw, m, T));
    else pf.reset(new PgenFile(path, n_raw, allele_ct, T));
    dosage = pf && pf->dosage();
  }
  // Bytes of the read buffers of one thread.
  double thread_bytes() const {
    const double rawb = static_cast<double>((static_cast<std::size_t>(n_raw) + 3) / 4);
    return (bf ? rawb : (dosage ? 3.0 * n_raw : rawb)) + 1048576.0;
  }
};

// The variants of a file for a scan of all of them: how many it has, the number of alleles of each (PGEN) and the
// file rows of those with two, which are the ones scanned. A BED says how many variants it has by its size (3 bytes,
// then (n + 3) / 4 bytes per variant) and all have two alleles; for a PGEN the .pvar is read for the ALT alleles of
// every line (variant_text.h).
struct FileVariants {
  std::size_t m = 0;
  std::vector<int> allele_ct, rows;
};
inline FileVariants file_variants(const std::string& path, bool bed, int n_raw, const std::string& variant_text) {
  FileVariants v;
  if (bed) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) Rcpp::stop("Cannot open BED file: " + path);
    const std::uint64_t size = static_cast<std::uint64_t>(f.tellg()), stride = (static_cast<std::uint64_t>(n_raw) + 3) / 4;
    if (size < 3 || stride == 0 || (size - 3) % stride != 0)
      Rcpp::stop("BED size does not match its FAM file: " + path);
    v.m = static_cast<std::size_t>((size - 3) / stride);
    v.rows.resize(v.m);
    for (std::size_t j = 0; j < v.m; ++j) v.rows[j] = static_cast<int>(j);
  } else {
    v.allele_ct = allele_counts(variant_text);
    v.m = v.allele_ct.size();
    for (std::size_t j = 0; j < v.m; ++j)
      if (v.allele_ct[j] == 2) v.rows.push_back(static_cast<int>(j));
  }
  if (v.rows.empty()) Rcpp::stop("No biallelic variant to test.");
  return v;
}

// The variants of one thread's block over the analysis samples, compact: 2-bit codes (hard calls) or 0.01 dosage steps,
// with the observed mean of each variant (a missing cell counts as that mean: 0 after centring) and whether it cannot be
// tested (monomorphic, all missing, or not read).
struct ScanBlock {
  int cap = 0, n = 0;
  bool dosage = false;
  std::size_t sb = 0;
  std::vector<unsigned char> g8;   // cap x sb: 2-bit codes (0/1/2 = count, 3 = missing)
  std::vector<std::uint16_t> g16;  // cap x n: dosage in steps of 0.01 (dosage_grid), PgenFile::kMissing = missing
  std::vector<double> mean;
  std::vector<int> mono;

  // Bytes of one variant of a block over n samples.
  static double var_bytes(int n, bool dos) {
    const double N = n;
    return dos ? 2.0 * N : (N + 3) / 4.0;
  }
  void init(int capacity, int n_, bool dos) {
    cap = capacity; n = n_; dosage = dos;
    sb = (static_cast<std::size_t>(n) + 3) / 4;
    if (dosage) g16.assign(static_cast<std::size_t>(cap) * n, 0); else g8.assign(cap * sb, 0);
    mean.assign(cap, 0.0); mono.assign(cap, 0);
  }
  // Centred values of rows row0 .. row0 + rows - 1 of slot s (row0 a multiple of 4).
  void centred(int s, int row0, int rows, double* out) const {
    const double mu = mean[s];
    if (dosage) {
      const std::uint16_t* d = &g16[static_cast<std::size_t>(s) * n + row0];
      for (int i = 0; i < rows; ++i) out[i] = d[i] == PgenFile::kMissing ? 0.0 : d[i] * 0.01 - mu;
      return;
    }
    const double lut[4] = {-mu, 1.0 - mu, 2.0 - mu, 0.0};
    const unsigned char* p = &g8[s * sb + (row0 >> 2)];
    int i = 0;
    for (; i + 4 <= rows; i += 4) {
      const unsigned b = p[i >> 2];
      out[i] = lut[b & 3]; out[i + 1] = lut[(b >> 2) & 3]; out[i + 2] = lut[(b >> 4) & 3]; out[i + 3] = lut[b >> 6];
    }
    for (; i < rows; ++i) out[i] = lut[(p[i >> 2] >> ((i & 3) << 1)) & 3];
  }
  // Centred value of sample i of slot s.
  double value(int s, int i) const {
    const double mu = mean[s];
    if (dosage) {
      const std::uint16_t v = g16[static_cast<std::size_t>(s) * n + i];
      return v == PgenFile::kMissing ? 0.0 : v * 0.01 - mu;
    }
    const unsigned code = (g8[s * sb + (i >> 2)] >> ((i & 3) << 1)) & 3;
    return code == 3 ? 0.0 : static_cast<double>(code) - mu;
  }
  // Reads variants j0 .. j0 + nb - 1 of the scan (thread th) into slots 0 .. nb - 1 and writes their observed count and
  // sum to n_obs[j] and a1_sum[j]. Returns the last variant whose read failed, or -1.
  long long read(ScanFile& F, long long j0, int nb, int th, double* n_obs, double* a1_sum) {
    long long bad = -1;
    for (int s = 0; s < nb; ++s) {
      const long long j = j0 + s;
      double no = 0, sum = 0, sq = 0;
      if (dosage) {
        const std::uint16_t* d = F.pf->read_dosage(static_cast<std::size_t>(F.vi[j]), th);
        if (!d) { bad = j; mono[s] = 1; continue; }
        std::uint16_t* o = &g16[static_cast<std::size_t>(s) * n];
        for (int i = 0; i < n; ++i) {
          const std::uint16_t v = d[F.si[i]];
          if (v == PgenFile::kMissing) { o[i] = v; continue; }
          const unsigned int c = dosage_grid(v);
          o[i] = static_cast<std::uint16_t>(c);
          const double x = c * 0.01;
          no += 1; sum += x; sq += x * x;
        }
      } else {
        const unsigned char* gv = F.bf ? F.bf->read_counts(static_cast<std::size_t>(F.vi[j]), th)
                                       : F.pf->read(static_cast<std::size_t>(F.vi[j]), th);
        if (!gv) { bad = j; mono[s] = 1; continue; }
        unsigned char* dst = &g8[s * sb];
        compact_genovec(gv, F.si.data(), n, F.ident, dst);
        long long cnt[4];
        count_codes(dst, sb, cnt);
        cnt[0] -= static_cast<long long>(sb * 4 - n);  // padding is code 0
        no = static_cast<double>(cnt[0] + cnt[1] + cnt[2]);
        sum = static_cast<double>(cnt[1] + 2 * cnt[2]);
        sq = static_cast<double>(cnt[1] + 4 * cnt[2]);
      }
      n_obs[j] = no; a1_sum[j] = sum;
      const double mu = no > 0 ? sum / no : 0.0, v = no > 0 ? sq / no - mu * mu : 0.0;
      mean[s] = mu;
      mono[s] = (no < 2 || !(v > 1e-12 * std::max(1.0, sq / std::max(no, 1.0)))) ? 1 : 0;
    }
    return bad;
  }
};

// Variants per thread: `block` when positive, else the most the memory budget allows (at most kScanBlockMax); never more
// than one pass over the p variants needs. shared = bytes held once, per_thread = bytes of a thread besides its variants,
// per_var = bytes per variant of a thread.
inline int scan_block_size(int block, double mem_gb, double shared, double per_thread, double per_var, int p, int ns,
                           int T) {
  const double GiB = 1073741824.0;
  const double budget = mem_gb * GiB * 0.95 - 268435456.0 - shared;
  const int per_thread_max = (p + T - 1) / T;
  int b = block;
  if (b <= 0) {
    const double room = budget / T - per_thread;
    if (room < per_var)
      Rcpp::stop("mem_gb = %.2f is too small: about %.2f GB are needed for %d samples and %d threads.", mem_gb,
                 (shared + T * (per_thread + per_var) + 268435456.0) / 0.95 / GiB, ns, T);
    b = static_cast<int>(std::min<double>(kScanBlockMax, std::floor(room / per_var)));
  }
  return std::max(1, std::min(b, std::max(1, per_thread_max)));
}

}  // namespace plink4r

#endif
