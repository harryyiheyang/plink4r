// PGEN reader on the bundled pgenlib (read side), one PgenReader per thread: packed 2-bit genovecs
// (read), for files whose header flags dosages 1/16384 fixed-point ALT dosages of every sample (read_dosage), and
// the file's own sparse form over a sample subset that pgenlib takes word by word (set_samples, read_sparse).
// Hard calls with their phase (read_phased): the genovec, the heterozygotes whose phase the file holds and which
// haplotype of each carries ALT.
// The one PgenFile of the package (geno objects, the association scans and the haplotypes).
// compact_genovec / count_codes: the sample subset of a genovec, and the calls of each code.
#ifndef PLINK4R_PGEN_READER_H
#define PLINK4R_PGEN_READER_H

#include <Rcpp.h>
#include <fstream>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string>
#include <cstdint>
#include <cstring>
#include <vector>
#include "pgenlib/pgenlib_read.h"

namespace plink4r {

// A file dosage (1/16384 steps, 0..32768) taken to the nearest 0.01: steps of 0.01, 0..200. PI 2026-11-18: the
// package holds dosages at 0.01, the precision imputation writes; 0.05 is 819/16384 in the file and gives 5.
inline unsigned int dosage_grid(std::uint16_t v) {
  return (static_cast<unsigned int>(v) * 100u + 8192u) >> 14;
}

// OpenMP thread number inside a parallel region (0 without OpenMP): the package's one copy.
inline int thread_id() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// Size of the current OpenMP team (1 without OpenMP; may be below num_threads(T) when the runtime limits it).
inline int thread_team() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

struct PgenHeader {
  plink2::PgenFileInfo info;
  unsigned char* alloc = nullptr;
  std::vector<uintptr_t> allele_idx_offsets, nonref_flags;
  PgenHeader() { plink2::PreinitPgfi(&info); }
  PgenHeader(const PgenHeader&) = delete;
  PgenHeader& operator=(const PgenHeader&) = delete;
  ~PgenHeader() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(&info, &e);
    if (alloc) plink2::aligned_free(alloc);
  }
};

struct ThreadReader {
  plink2::PgenReader pgr;
  plink2::PgrSampleSubsetIndex pssi;
  unsigned char* alloc = nullptr;
  unsigned char* genovec = nullptr;
  unsigned char* dosage_present = nullptr;
  unsigned char* dosage_main = nullptr;
  unsigned char* phasepresent = nullptr;
  unsigned char* phaseinfo = nullptr;
  std::vector<std::uint16_t> full;
  ThreadReader() { plink2::PreinitPgr(&pgr); }
  ThreadReader(const ThreadReader&) = delete;
  ThreadReader& operator=(const ThreadReader&) = delete;
  ~ThreadReader() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgr(&pgr, &e);
    if (alloc) plink2::aligned_free(alloc);
    if (genovec) plink2::aligned_free(genovec);
    if (dosage_present) plink2::aligned_free(dosage_present);
    if (dosage_main) plink2::aligned_free(dosage_main);
    if (phasepresent) plink2::aligned_free(phasepresent);
    if (phaseinfo) plink2::aligned_free(phaseinfo);
  }
};

class PgenFile {
public:
  PgenFile(const std::string& path, std::size_t n, const std::vector<int>& allele_ct, int threads)
    : path_(path) {
    plink2::PgenFileInfo& info = header_.info;
    char err[plink2::kPglErrstrBufBlen];
    err[0] = '\0';
    plink2::PgenHeaderCtrl header_ctrl;
    uintptr_t pgfi_cachelines;
    if (plink2::PgfiInitPhase1(path.c_str(), nullptr, UINT32_MAX, UINT32_MAX, &header_ctrl,
                               &info, &pgfi_cachelines, err) != plink2::kPglRetSuccess)
      fail(err);
    if (info.raw_sample_ct != n) Rcpp::stop("PGEN sample count does not match .psam: " + path);
    const std::size_t m = info.raw_variant_ct;
    if (allele_ct.size() != m) Rcpp::stop("PGEN variant count does not match .pvar: " + path);
    uint32_t max_allele_ct = 2;
    for (int a : allele_ct) max_allele_ct = std::max<uint32_t>(max_allele_ct, a);
    info.max_allele_ct = 2;
    if (max_allele_ct > 2 || (header_ctrl & 0x30)) {
      std::vector<uintptr_t>& offsets = header_.allele_idx_offsets;
      offsets.resize(m + 1);
      offsets[0] = 0;
      for (std::size_t j = 0; j < m; ++j) offsets[j + 1] = offsets[j] + allele_ct[j];
      info.allele_idx_offsets = offsets.data();
      info.max_allele_ct = max_allele_ct;
    }
    if ((header_ctrl & 0xc0) == 0xc0) {
      header_.nonref_flags.resize(plink2::DivUp(m, plink2::kBitsPerWord) + 1);
      info.nonref_flags = header_.nonref_flags.data();
    }
    if (plink2::cachealigned_malloc(std::max<uintptr_t>(pgfi_cachelines, 1) * plink2::kCacheline, &header_.alloc))
      Rcpp::stop("Out of memory opening PGEN: " + path);
    uint32_t max_vrec_width;
    uintptr_t pgr_cachelines;
    if (plink2::PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, info.raw_variant_ct, &max_vrec_width,
                               &info, header_.alloc, &pgr_cachelines, err) != plink2::kPglRetSuccess)
      fail(err);
    genovec_bytes_ = plink2::NypCtToVecCt(n) * plink2::kBytesPerVec;
    dosage_ = (info.gflags & plink2::kfPgenGlobalDosagePresent) != 0;
    const uintptr_t present_bytes = ((n + 7) / 8 / plink2::kBytesPerVec + 2) * plink2::kBytesPerVec;
    const uintptr_t main_bytes = (2 * n / plink2::kBytesPerVec + 2) * plink2::kBytesPerVec;
    for (int t = 0; t < threads; ++t) {
      readers_.emplace_back(new ThreadReader);
      ThreadReader& r = *readers_.back();
      plink2::PgrSetFreadBuf(nullptr, &r.pgr);
      if (plink2::cachealigned_malloc(std::max<uintptr_t>(pgr_cachelines, 1) * plink2::kCacheline, &r.alloc) ||
          plink2::cachealigned_malloc(genovec_bytes_, &r.genovec) ||
          (dosage_ && (plink2::cachealigned_malloc(present_bytes, &r.dosage_present) ||
                       plink2::cachealigned_malloc(main_bytes, &r.dosage_main))))
        Rcpp::stop("Out of memory opening PGEN: " + path);
      if (dosage_) r.full.resize(n);
      if (plink2::PgrInit(path.c_str(), max_vrec_width, &info, &r.pgr, r.alloc) != plink2::kPglRetSuccess)
        Rcpp::stop("Could not open PGEN reader: " + path);
      plink2::PgrClearSampleSubsetIndex(&r.pgr, &r.pssi);
    }
  }
  std::size_t samples() const { return header_.info.raw_sample_ct; }
  // Packed 2-bit genovec (all samples, file order); codes 0/1/2 = ALT count, 3 = missing.
  const unsigned char* read(std::size_t variant, int thread) {
    ThreadReader& r = *readers_[thread];
    if (plink2::PgrGet(nullptr, r.pssi, static_cast<uint32_t>(samples()), static_cast<uint32_t>(variant),
                       &r.pgr, reinterpret_cast<uintptr_t*>(r.genovec)) != plink2::kPglRetSuccess)
      return nullptr;
    return r.genovec;
  }

  bool dosage() const { return dosage_; }

  // One variant with its phase, all samples in file order: genovec (codes 0/1/2 = ALT count, 3 = missing; a
  // multiallelic variant counts every non-REF allele), present = the bit of each heterozygote whose phase the
  // file holds (present_ct of them), info = for those, bit set: ALT on the first haplotype (1|0), clear: on the
  // second (0|1). The bits of info outside present mean nothing. False on a read error.
  struct Phased {
    const std::uint64_t* genovec = nullptr;
    const std::uint64_t* present = nullptr;
    const std::uint64_t* info = nullptr;
    uint32_t present_ct = 0;
  };
  bool read_phased(std::size_t variant, int thread, Phased& out) {
    ThreadReader& r = *readers_[thread];
    if (!r.phasepresent) {
      const uintptr_t bytes = plink2::BitCtToVecCt(samples()) * plink2::kBytesPerVec;
      if (plink2::cachealigned_malloc(bytes, &r.phasepresent) || plink2::cachealigned_malloc(bytes, &r.phaseinfo))
        return false;
    }
    out.present_ct = 0;
    if (plink2::PgrGetP(nullptr, r.pssi, static_cast<uint32_t>(samples()), static_cast<uint32_t>(variant), &r.pgr,
                        reinterpret_cast<uintptr_t*>(r.genovec), reinterpret_cast<uintptr_t*>(r.phasepresent),
                        reinterpret_cast<uintptr_t*>(r.phaseinfo), &out.present_ct) != plink2::kPglRetSuccess)
      return false;
    out.genovec = reinterpret_cast<const std::uint64_t*>(r.genovec);
    out.present = reinterpret_cast<const std::uint64_t*>(r.phasepresent);
    out.info = reinterpret_cast<const std::uint64_t*>(r.phaseinfo);
    return true;
  }

  // The samples read_sparse() returns: rows (0-based .psam rows, ascending, distinct); empty or every row keeps
  // all samples. pgenlib subsets the genovec and the dosages word by word.
  void set_samples(const std::vector<int>& rows) {
    const std::size_t n = samples();
    include_.clear();
    cumpop_.clear();
    sample_ct_ = static_cast<uint32_t>(rows.empty() ? n : rows.size());
    if (rows.empty() || rows.size() == n) {
      for (auto& r : readers_) plink2::PgrClearSampleSubsetIndex(&r->pgr, &r->pssi);
      return;
    }
    const uint32_t wc = static_cast<uint32_t>(plink2::BitCtToWordCt(n));
    include_.assign(wc, 0);
    for (int i : rows) include_[i / plink2::kBitsPerWord] |= uintptr_t(1) << (i % plink2::kBitsPerWord);
    cumpop_.assign(wc, 0);
    plink2::FillCumulativePopcounts(include_.data(), wc, cumpop_.data());
    for (auto& r : readers_) plink2::PgrSetSampleSubsetIndex(cumpop_.data(), &r->pgr, &r->pssi);
  }

  // One variant over the samples of set_samples() in the file's sparse form: genovec (codes 0/1/2 = ALT count,
  // 3 = missing, 32 samples per word) and, for a dosage file, the bit of each sample with an explicit dosage
  // (present) and those dosages in sample order (main, 1/16384 steps; dosage_ct of them). An explicit dosage
  // overrides the sample's genovec code. False on a read error.
  struct Sparse {
    const std::uint64_t* genovec = nullptr;
    const std::uint64_t* present = nullptr;
    const std::uint16_t* main = nullptr;
    uint32_t dosage_ct = 0;
  };
  bool read_sparse(std::size_t variant, int thread, Sparse& out) {
    ThreadReader& r = *readers_[thread];
    const uintptr_t* inc = include_.empty() ? nullptr : include_.data();
    const uint32_t n = sample_ct_ ? sample_ct_ : static_cast<uint32_t>(samples());
    out.dosage_ct = 0;
    if (dosage_) {
      if (plink2::PgrGetD(inc, r.pssi, n, static_cast<uint32_t>(variant), &r.pgr,
                          reinterpret_cast<uintptr_t*>(r.genovec), reinterpret_cast<uintptr_t*>(r.dosage_present),
                          reinterpret_cast<std::uint16_t*>(r.dosage_main), &out.dosage_ct) != plink2::kPglRetSuccess)
        return false;
    } else if (plink2::PgrGet(inc, r.pssi, n, static_cast<uint32_t>(variant), &r.pgr,
                              reinterpret_cast<uintptr_t*>(r.genovec)) != plink2::kPglRetSuccess) {
      return false;
    }
    out.genovec = reinterpret_cast<const std::uint64_t*>(r.genovec);
    out.present = reinterpret_cast<const std::uint64_t*>(r.dosage_present);
    out.main = reinterpret_cast<const std::uint16_t*>(r.dosage_main);
    return true;
  }
  // ALT dosage of every sample (file order) in 1/16384 units (0..32768), kMissing for a missing call;
  // explicit dosages override hardcalls (adapted from SuSiE4I pgen_file.h).
  static constexpr std::uint16_t kMissing = 65535;
  const std::uint16_t* read_dosage(std::size_t variant, int thread) {
    ThreadReader& r = *readers_[thread];
    const uint32_t n = static_cast<uint32_t>(samples());
    uint32_t dosage_ct = 0;
    if (plink2::PgrGetD(nullptr, r.pssi, n, static_cast<uint32_t>(variant), &r.pgr,
                        reinterpret_cast<uintptr_t*>(r.genovec), reinterpret_cast<uintptr_t*>(r.dosage_present),
                        reinterpret_cast<std::uint16_t*>(r.dosage_main), &dosage_ct) != plink2::kPglRetSuccess)
      return nullptr;
    static const std::uint16_t code_dosage[4] = {0, 16384, 32768, kMissing};
    std::uint16_t* full = r.full.data();
    for (uint32_t i = 0; i < n; ++i) full[i] = code_dosage[(r.genovec[i >> 2] >> (2 * (i & 3u))) & 3u];
    if (dosage_ct) {
      const std::uint64_t* present = reinterpret_cast<const std::uint64_t*>(r.dosage_present);
      const std::uint16_t* main = reinterpret_cast<const std::uint16_t*>(r.dosage_main);
      uint32_t k = 0;
      for (uint32_t w = 0; w < (n + 63) / 64 && k < dosage_ct; ++w)
        for (std::uint64_t bits = present[w]; bits && k < dosage_ct; bits &= bits - 1)
          full[64 * w + __builtin_ctzll(bits)] = main[k++];
    }
    return full;
  }

private:
  [[noreturn]] void fail(const char* err) const {
    std::string msg(err);
    if (msg.compare(0, 7, "Error: ") == 0) msg.erase(0, 7);
    while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) msg.pop_back();
    Rcpp::stop("Could not open PGEN file " + path_ + (msg.empty() ? "" : ": " + msg));
  }
  bool dosage_ = false;
  uint32_t sample_ct_ = 0;
  std::vector<uintptr_t> include_;
  std::vector<uint32_t> cumpop_;
  std::string path_;
  uintptr_t genovec_bytes_ = 0;
  PgenHeader header_;
  std::vector<std::unique_ptr<ThreadReader>> readers_;
};

// Genovec of the samples sidx[0..ns-1] (0-based .psam rows, output order) into dst ((ns + 3) / 4 bytes): codes
// 0/1/2/3 as in the file; the bits past the last sample are code 0, so they add nothing to counts or products.
// identity: sidx is 0..ns-1 over all n_raw samples, a plain copy.
inline void compact_genovec(const unsigned char* gv, const int* sidx, int ns, bool identity, unsigned char* dst) {
  const std::size_t sb = (static_cast<std::size_t>(ns) + 3) / 4;
  if (identity) {
    std::memcpy(dst, gv, sb);
    const int pad = static_cast<int>(sb * 4 - ns);
    if (pad) dst[sb - 1] &= static_cast<unsigned char>(0xFFu >> (2 * pad));
    return;
  }
  std::memset(dst, 0, sb);
  for (int s = 0; s < ns; ++s) {
    const int r = sidx[s];
    dst[s >> 2] |= static_cast<unsigned char>(((gv[r >> 2] >> (2 * (r & 3))) & 3u) << (2 * (s & 3)));
  }
}

// True when sidx is 0..n_raw-1 (every sample, file order): compact_genovec is then a copy.
inline bool identity_samples(const int* sidx, int ns, int n_raw) {
  if (ns != n_raw) return false;
  for (int s = 0; s < ns; ++s) if (sidx[s] != s) return false;
  return true;
}

// Codes 0/1/2/3 of a compacted genovec of sb bytes (padding is code 0): cnt[c] = calls with code c.
inline void count_codes(const unsigned char* g, std::size_t sb, long long cnt[4]) {
  cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
  for (std::size_t y = 0; y < sb; ++y) {
    const unsigned v = g[y];
    ++cnt[v & 3u]; ++cnt[(v >> 2) & 3u]; ++cnt[(v >> 4) & 3u]; ++cnt[v >> 6];
  }
}

// PLINK 1 SNP-major .bed: header and size checked on opening, one stream per thread. read() returns the row of
// a variant (all samples, file order) in the .bed coding (00 two A1, 01 missing, 10 one A1, 11 no A1);
// read_counts() the same row recoded as a genovec counting A1 (codes 0/1/2, 3 missing), so the PGEN helpers
// (compact_genovec, count_codes) apply. A stream is moved only when the variant asked for is not the one after the
// last it read: a seek empties the buffer of the stream, so a seek per variant costs a read of a whole buffer for
// a row that may be a few bytes (10 million rows of 5 bytes: 4.8 s with a seek each, 0.13 s without). The package's
// one .bed reader.
class BedFile {
public:
  BedFile(const std::string& path, std::size_t n, std::size_t m, int threads)
    : path_(path), stride_((n + 3) / 4) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) Rcpp::stop("Cannot open BED file: " + path);
    const std::streamoff size = f.tellg();
    if (m > 0 && (size < 0 || static_cast<std::uint64_t>(size) != 3 + static_cast<std::uint64_t>(m) * stride_))
      Rcpp::stop("BED size does not match its BIM and FAM files: " + path);
    f.seekg(0);
    unsigned char h[3];
    if (!f.read(reinterpret_cast<char*>(h), 3) || h[0] != 0x6c || h[1] != 0x1b || h[2] != 0x01)
      Rcpp::stop("BED must be a PLINK SNP-major BED file: " + path);
    for (int t = 0; t < std::max(1, threads); ++t) {
      readers_.emplace_back(new Reader);
      readers_.back()->f.open(path, std::ios::binary);
      readers_.back()->row.resize(stride_);
    }
  }
  std::size_t stride() const { return stride_; }
  const unsigned char* read(std::size_t variant, int thread) { return load(variant, thread); }
  const unsigned char* read_counts(std::size_t variant, int thread) {
    unsigned char* r = load(variant, thread);
    if (!r) return nullptr;
    to_counts(r, stride_);
    return r;
  }
  // A .bed row recoded in place as a genovec counting A1 (codes 0/1/2, 3 missing).
  static void to_counts(unsigned char* row, std::size_t stride) {
    struct Tab { unsigned char t[256];
      Tab() { static const unsigned char m[4] = {2, 3, 1, 0};
        for (int b = 0; b < 256; ++b) { unsigned v = 0;
          for (int k = 0; k < 4; ++k) v |= static_cast<unsigned>(m[(b >> (2 * k)) & 3]) << (2 * k);
          t[b] = static_cast<unsigned char>(v); } } };
    static const Tab T;
    for (std::size_t y = 0; y < stride; ++y) row[y] = T.t[row[y]];
  }
private:
  // The stream of one thread, the row it read last and the place in the file where it stands (kNowhere: not known).
  static constexpr std::uint64_t kNowhere = ~static_cast<std::uint64_t>(0);
  struct Reader {
    std::ifstream f;
    std::vector<unsigned char> row;
    std::uint64_t at = kNowhere;
  };
  unsigned char* load(std::size_t variant, int thread) {
    Reader& r = *readers_[thread];
    const std::uint64_t pos = 3 + static_cast<std::uint64_t>(variant) * stride_;
    if (r.at != pos) r.f.seekg(static_cast<std::streamoff>(pos));
    if (!r.f.read(reinterpret_cast<char*>(r.row.data()), static_cast<std::streamsize>(stride_))) {
      r.f.clear();
      r.at = kNowhere;
      return nullptr;
    }
    r.at = pos + stride_;
    return r.row.data();
  }
  std::string path_;
  std::size_t stride_;
  std::vector<std::unique_ptr<Reader>> readers_;
};

}  // namespace plink4r
#endif
