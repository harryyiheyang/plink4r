// PGEN reader on the bundled pgenlib (read side), one PgenReader per thread. A variant is read as (file line, ALT):
// the count of the ALT-th ALT allele of that line. A line with one ALT allele is read as pgenlib reads a biallelic
// variant (PgrGet, PgrGetD, PgrGetP); a line with several is read by allele (PgrGet1, PgrGet1D), a sample that
// carries another ALT of the line counting 0, and its two haplotypes from the allele codes of the line (PgrGetMP).
// What is read: packed 2-bit genovecs (read), for files whose header flags dosages 1/16384 fixed-point dosages of
// every sample (read_dosage), the file's own sparse form over a sample subset that pgenlib takes word by word
// (set_samples, read_sparse), and the allele on each of the two haplotypes of every sample (read_haplotypes).
// The one PgenFile of the package (geno objects, the association scans and the haplotypes).
// compact_genovec / count_codes: the sample subset of a genovec, and the calls of each code.
#ifndef PLINK4R_PGEN_READER_H
#define PLINK4R_PGEN_READER_H

#include <Rcpp.h>
#include <atomic>
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

namespace plink2 {
// pgenlib_read.cc: the cachelines PgrInit() takes from a reader's block for a file with these flags. PgfiInitPhase2()
// reports fewer for a file of fixed-width records (storage modes 0x02 to 0x04): only the base genovec and the read
// buffer, while PgrInit() also places the dosage and phase workspaces there.
uintptr_t CountPgrAllocCachelinesRequired(uint32_t raw_sample_ct, PgenGlobalFlags gflags, uint32_t max_allele_ct,
                                          uint32_t fread_buf_byte_ct);
}

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

// A block for pgenlib (cacheline-aligned), freed with its owner.
struct AlignedBlock {
  unsigned char* p = nullptr;
  AlignedBlock() = default;
  AlignedBlock(const AlignedBlock&) = delete;
  AlignedBlock& operator=(const AlignedBlock&) = delete;
  ~AlignedBlock() { if (p) plink2::aligned_free(p); }
  bool take(uintptr_t bytes) { return !plink2::cachealigned_malloc(std::max<uintptr_t>(bytes, 1), &p); }
  uintptr_t* words() const { return reinterpret_cast<uintptr_t*>(p); }
};

struct ThreadReader {
  plink2::PgenReader pgr;
  plink2::PgrSampleSubsetIndex pssi;
  AlignedBlock alloc, genovec, dosage_present, dosage_main;
  std::vector<std::uint16_t> full;
  // haplotypes (read_haplotypes): phase and the allele codes of a line with several ALT alleles as pgenlib gives
  // them, and the line they were last decoded for into h1, h2 (allele on each haplotype), state
  AlignedBlock phasepresent, phaseinfo, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals;
  std::vector<unsigned char> h1, h2, state;
  long long hap_line = -1;
  ThreadReader() { plink2::PreinitPgr(&pgr); }
  ThreadReader(const ThreadReader&) = delete;
  ThreadReader& operator=(const ThreadReader&) = delete;
  ~ThreadReader() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgr(&pgr, &e);
  }
};

class PgenFile {
public:
  // allele_ct: the number of alleles of every file line, as its .pvar names them.
  PgenFile(const std::string& path, std::size_t n, const std::vector<int>& allele_ct, int threads)
    : path_(path), allele_ct_(allele_ct) {
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
    pgr_cachelines = std::max(pgr_cachelines, plink2::CountPgrAllocCachelinesRequired(
        info.raw_sample_ct, info.gflags, info.max_allele_ct, max_vrec_width));
    genovec_bytes_ = plink2::NypCtToVecCt(n) * plink2::kBytesPerVec;
    dosage_ = (info.gflags & plink2::kfPgenGlobalDosagePresent) != 0;
    const uintptr_t present_bytes = ((n + 7) / 8 / plink2::kBytesPerVec + 2) * plink2::kBytesPerVec;
    const uintptr_t main_bytes = (2 * n / plink2::kBytesPerVec + 2) * plink2::kBytesPerVec;
    for (int t = 0; t < threads; ++t) {
      readers_.emplace_back(new ThreadReader);
      ThreadReader& r = *readers_.back();
      plink2::PgrSetFreadBuf(nullptr, &r.pgr);
      if (!r.alloc.take(std::max<uintptr_t>(pgr_cachelines, 1) * plink2::kCacheline) || !r.genovec.take(genovec_bytes_) ||
          (dosage_ && (!r.dosage_present.take(present_bytes) || !r.dosage_main.take(main_bytes))))
        Rcpp::stop("Out of memory opening PGEN: " + path);
      if (dosage_) r.full.resize(n);
      if (plink2::PgrInit(path.c_str(), max_vrec_width, &info, &r.pgr, r.alloc.p) != plink2::kPglRetSuccess)
        Rcpp::stop("Could not open PGEN reader: " + path);
      plink2::PgrClearSampleSubsetIndex(&r.pgr, &r.pssi);
    }
  }
  std::size_t samples() const { return header_.info.raw_sample_ct; }
  // Whether file line `line` names several ALT alleles.
  bool several(std::size_t line) const { return allele_ct_[line] > 2; }

  // Packed 2-bit genovec of (line, alt) (all samples, file order); codes 0/1/2 = count of that ALT, 3 = missing.
  const unsigned char* read(std::size_t line, int alt, int thread) {
    ThreadReader& r = *readers_[thread];
    const uint32_t n = static_cast<uint32_t>(samples()), v = static_cast<uint32_t>(line);
    const plink2::PglErr e = several(line)
        ? plink2::PgrGet1(nullptr, r.pssi, n, v, static_cast<uint32_t>(alt), &r.pgr, r.genovec.words())
        : plink2::PgrGet(nullptr, r.pssi, n, v, &r.pgr, r.genovec.words());
    return e == plink2::kPglRetSuccess ? r.genovec.p : nullptr;
  }

  bool dosage() const { return dosage_; }

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

  // (line, alt) over the samples of set_samples() in the file's sparse form: genovec (codes 0/1/2 = count of that
  // ALT, 3 = missing, 32 samples per word) and, for a dosage file, the bit of each sample with an explicit dosage
  // (present) and those dosages in sample order (main, 1/16384 steps; dosage_ct of them). An explicit dosage
  // overrides the sample's genovec code. False on a read error (failure() says which).
  struct Sparse {
    const std::uint64_t* genovec = nullptr;
    const std::uint64_t* present = nullptr;
    const std::uint16_t* main = nullptr;
    uint32_t dosage_ct = 0;
  };
  bool read_sparse(std::size_t line, int alt, int thread, Sparse& out) {
    ThreadReader& r = *readers_[thread];
    const uintptr_t* inc = include_.empty() ? nullptr : include_.data();
    const uint32_t n = sample_ct_ ? sample_ct_ : static_cast<uint32_t>(samples());
    out.dosage_ct = 0;
    if (!get(r, inc, n, line, alt, &out.dosage_ct)) return false;
    out.genovec = reinterpret_cast<const std::uint64_t*>(r.genovec.p);
    out.present = reinterpret_cast<const std::uint64_t*>(r.dosage_present.p);
    out.main = reinterpret_cast<const std::uint16_t*>(r.dosage_main.p);
    return true;
  }
  // Dosage of the ALT of (line, alt) of every sample (file order) in 1/16384 units (0..32768), kMissing for a
  // missing call; explicit dosages override hardcalls (adapted from SuSiE4I pgen_file.h).
  static constexpr std::uint16_t kMissing = 65535;
  const std::uint16_t* read_dosage(std::size_t line, int alt, int thread) {
    ThreadReader& r = *readers_[thread];
    const uint32_t n = static_cast<uint32_t>(samples());
    uint32_t dosage_ct = 0;
    if (!get(r, nullptr, n, line, alt, &dosage_ct)) return nullptr;
    static const std::uint16_t code_dosage[4] = {0, 16384, 32768, kMissing};
    std::uint16_t* full = r.full.data();
    const unsigned char* gv = r.genovec.p;
    for (uint32_t i = 0; i < n; ++i) full[i] = code_dosage[(gv[i >> 2] >> (2 * (i & 3u))) & 3u];
    if (dosage_ct) {
      const std::uint64_t* present = reinterpret_cast<const std::uint64_t*>(r.dosage_present.p);
      const std::uint16_t* main = reinterpret_cast<const std::uint16_t*>(r.dosage_main.p);
      uint32_t k = 0;
      for (uint32_t w = 0; w < (n + 63) / 64 && k < dosage_ct; ++w)
        for (std::uint64_t bits = present[w]; bits && k < dosage_ct; bits &= bits - 1)
          full[64 * w + __builtin_ctzll(bits)] = main[k++];
    }
    return full;
  }

  // The allele on each of the two haplotypes of every sample at file line `line` (all samples, file order): h1 and
  // h2 hold allele codes (0 = REF, k = the k-th ALT of the line; a homozygote the same code twice), state says what
  // is known: kPhased (h1 is the first haplotype, h2 the second, as the file phases them; every homozygote),
  // kUnphased (a heterozygote whose phase the file does not hold: h1 < h2 are its two alleles, on which haplotype
  // each is unknown) or kNoCall (missing call). A thread keeps the line it decoded last, so the ALT alleles of one
  // line are served by one read. False on a read error.
  enum : unsigned char { kPhased = 0, kUnphased = 1, kNoCall = 2 };
  struct Haplotypes {
    const unsigned char* h1 = nullptr;
    const unsigned char* h2 = nullptr;
    const unsigned char* state = nullptr;
  };
  bool read_haplotypes(std::size_t line, int thread, Haplotypes& out) {
    ThreadReader& r = *readers_[thread];
    const uint32_t n = static_cast<uint32_t>(samples());
    out.h1 = r.h1.data(); out.h2 = r.h2.data(); out.state = r.state.data();
    if (r.hap_line == static_cast<long long>(line)) return true;
    if (r.h1.empty()) {
      const uintptr_t bits = plink2::BitCtToVecCt(n) * plink2::kBytesPerVec;
      const uintptr_t codes = (2 * static_cast<uintptr_t>(n) / plink2::kBytesPerVec + 2) * plink2::kBytesPerVec;
      if (!r.phasepresent.take(bits) || !r.phaseinfo.take(bits) || !r.patch_01_set.take(bits) ||
          !r.patch_10_set.take(bits) || !r.patch_01_vals.take(codes) || !r.patch_10_vals.take(codes))
        return false;
      r.h1.resize(n); r.h2.resize(n); r.state.resize(n);
      out.h1 = r.h1.data(); out.h2 = r.h2.data(); out.state = r.state.data();
    }
    r.hap_line = -1;
    plink2::PgenVariant pv;
    plink2::PreinitPgv(&pv);
    pv.genovec = r.genovec.words();
    pv.phasepresent = r.phasepresent.words(); pv.phaseinfo = r.phaseinfo.words();
    pv.patch_01_set = r.patch_01_set.words(); pv.patch_01_vals = r.patch_01_vals.p;
    pv.patch_10_set = r.patch_10_set.words(); pv.patch_10_vals = r.patch_10_vals.p;
    pv.patch_01_ct = pv.patch_10_ct = pv.phasepresent_ct = 0;
    const uint32_t v = static_cast<uint32_t>(line);
    const plink2::PglErr e = several(line)
        ? plink2::PgrGetMP(nullptr, r.pssi, n, v, &r.pgr, &pv)
        : plink2::PgrGetP(nullptr, r.pssi, n, v, &r.pgr, pv.genovec, pv.phasepresent, pv.phaseinfo, &pv.phasepresent_ct);
    if (e != plink2::kPglRetSuccess) return false;
    const std::uint64_t* present = reinterpret_cast<const std::uint64_t*>(r.phasepresent.p);
    const std::uint64_t* info = reinterpret_cast<const std::uint64_t*>(r.phaseinfo.p);
    const std::uint64_t* set01 = reinterpret_cast<const std::uint64_t*>(r.patch_01_set.p);
    const std::uint64_t* set10 = reinterpret_cast<const std::uint64_t*>(r.patch_10_set.p);
    const unsigned char* gv = r.genovec.p;
    unsigned char* h1 = r.h1.data(); unsigned char* h2 = r.h2.data(); unsigned char* st = r.state.data();
    uint32_t k01 = 0, k10 = 0;  // the patches come in sample order
    for (uint32_t i = 0; i < n; ++i) {
      const unsigned code = (gv[i >> 2] >> (2 * (i & 3u))) & 3u;
      if (code == 3) { h1[i] = h2[i] = 0; st[i] = kNoCall; continue; }
      unsigned char lo = 0, hi = 0;  // the two alleles, lo <= hi
      if (code == 1) {
        hi = 1;
        if (pv.patch_01_ct && ((set01[i >> 6] >> (i & 63u)) & 1u)) hi = pv.patch_01_vals[k01++];
      } else if (code == 2) {
        lo = hi = 1;
        if (pv.patch_10_ct && ((set10[i >> 6] >> (i & 63u)) & 1u)) {
          lo = pv.patch_10_vals[2 * k10]; hi = pv.patch_10_vals[2 * k10 + 1]; ++k10;
        }
      }
      if (lo == hi) { h1[i] = h2[i] = lo; st[i] = kPhased; continue; }
      if (pv.phasepresent_ct && ((present[i >> 6] >> (i & 63u)) & 1u)) {
        // the phase bit of a heterozygote: set = the allele of the higher code is on the first haplotype (1|0)
        const bool first_hi = ((info[i >> 6] >> (i & 63u)) & 1u) != 0;
        h1[i] = first_hi ? hi : lo; h2[i] = first_hi ? lo : hi; st[i] = kPhased;
      } else { h1[i] = lo; h2[i] = hi; st[i] = kUnphased; }
    }
    r.hap_line = static_cast<long long>(line);
    return true;
  }

  // What to say of a read that failed at file line `line` (0-based).
  std::string failure(std::size_t line) const {
    if (unsupported_.load())
      return "Line " + std::to_string(line + 1) + " of the variants of " + path_ + " names several ALT alleles and " +
             "its record holds dosages: pgenlib does not read such a record allele by allele.";
    return "Genotype read failed at line " + std::to_string(line + 1) + " of the variants of " + path_;
  }

private:
  // The hard calls of (line, alt) into the genovec of r and, for a dosage file, its explicit dosages into the
  // dosage buffers of r (*dosage_ct of them); inc, n: the sample subset (null: all n samples).
  bool get(ThreadReader& r, const uintptr_t* inc, uint32_t n, std::size_t line, int alt, uint32_t* dosage_ct) {
    const uint32_t v = static_cast<uint32_t>(line);
    plink2::PglErr e;
    if (!several(line)) {
      e = dosage_ ? plink2::PgrGetD(inc, r.pssi, n, v, &r.pgr, r.genovec.words(), r.dosage_present.words(),
                                    reinterpret_cast<std::uint16_t*>(r.dosage_main.p), dosage_ct)
                  : plink2::PgrGet(inc, r.pssi, n, v, &r.pgr, r.genovec.words());
    } else if (dosage_) {
      e = plink2::PgrGet1D(inc, r.pssi, n, v, static_cast<plink2::AlleleCode>(alt), &r.pgr, r.genovec.words(),
                           r.dosage_present.words(), reinterpret_cast<std::uint16_t*>(r.dosage_main.p), dosage_ct);
      if (e == plink2::kPglRetNotYetSupported) unsupported_.store(true);
    } else {
      e = plink2::PgrGet1(inc, r.pssi, n, v, static_cast<uint32_t>(alt), &r.pgr, r.genovec.words());
    }
    return e == plink2::kPglRetSuccess;
  }
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
  std::vector<int> allele_ct_;
  std::atomic<bool> unsupported_{false};
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
