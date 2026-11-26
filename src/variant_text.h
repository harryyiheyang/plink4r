// The variant text of a genotype file, read in compiled code: the .bim of a BED or the .pvar of a PGEN, line by line
// in file order. R holds none of it. Three layers:
//   TextLines      the lines of a text file, one after another
//   VariantLines   the variant lines of a .bim / .pvar, and what a field of such a line says (chr_code, parse_bp,
//                  allele_count, variant_fields)
//   what a scan takes from them: allele_counts (a PGEN needs the number of alleles of every variant before any
//   genotype is read), read_variant_table (chromosome, position and alleles of every line, for a scan that places
//   windows) and ResultText (the result of a scan of every variant: each line as it is, then its statistics).
#ifndef PLINK4R_VARIANT_TEXT_H
#define PLINK4R_VARIANT_TEXT_H

#include <Rcpp.h>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

namespace plink4r {

// The lines of a text file in order, without their end of line ("\n" or "\r\n").
class TextLines {
 public:
  explicit TextLines(const std::string& path) : path_(path), buf_(std::size_t(1) << 22) {
    f_ = std::fopen(path.c_str(), "rb");
    if (!f_) Rcpp::stop("Cannot open " + path);
  }
  ~TextLines() { if (f_) std::fclose(f_); }
  TextLines(const TextLines&) = delete;
  TextLines& operator=(const TextLines&) = delete;
  const std::string& path() const { return path_; }
  // The next line, s .. s + len - 1, valid until the next call; false at the end of the file.
  bool next(const char*& s, std::size_t& len) {
    for (;;) {
      const char* b = buf_.data() + pos_;
      const char* e = static_cast<const char*>(std::memchr(b, '\n', end_ - pos_));
      if (e) { s = b; len = static_cast<std::size_t>(e - b); pos_ += len + 1; break; }
      if (eof_) {  // the last line has no "\n"
        if (pos_ == end_) return false;
        s = b; len = end_ - pos_; pos_ = end_; break;
      }
      refill();
    }
    if (len && s[len - 1] == '\r') --len;
    return true;
  }

 private:
  void refill() {
    const std::size_t keep = end_ - pos_;
    if (pos_ > 0) std::memmove(buf_.data(), buf_.data() + pos_, keep);
    else if (keep == buf_.size()) buf_.resize(buf_.size() * 2);  // a line longer than the buffer
    pos_ = 0; end_ = keep;
    const std::size_t got = std::fread(buf_.data() + end_, 1, buf_.size() - end_, f_);
    if (got == 0) {
      if (std::ferror(f_)) Rcpp::stop("Read error in " + path_);
      eof_ = true;
    }
    end_ += got;
  }
  std::string path_;
  std::vector<char> buf_;
  std::FILE* f_ = nullptr;
  std::size_t pos_ = 0, end_ = 0;
  bool eof_ = false;
};

// The variant lines of a .bim (six columns, no header) or of a .pvar (its "##" lines are passed over and its
// "#CHROM" line is kept as header()).
class VariantLines {
 public:
  VariantLines(const std::string& path, bool bim) : in_(path) {
    if (bim) return;
    const char* s; std::size_t len;
    while (in_.next(s, len)) {
      if (len >= 2 && s[0] == '#' && s[1] == '#') continue;
      if (len >= 6 && std::memcmp(s, "#CHROM", 6) == 0) { header_.assign(s, len); return; }
      break;
    }
    Rcpp::stop("No #CHROM line before the variants of " + path);
  }
  const std::string& header() const { return header_; }
  const std::string& path() const { return in_.path(); }
  bool next(const char*& s, std::size_t& len) { return in_.next(s, len); }

 private:
  TextLines in_;
  std::string header_;
};

// The fields of a variant line that a scan reads, by their place (fields are separated by tabs or blanks):
// CHROM POS . REF ALT of a .pvar; chromosome . . position A1 A2 of a .bim, of which A1 (column 5) is the counted
// allele (alt) and A2 the other (ref).
struct VariantFields {
  const char* chr = nullptr; const char* pos = nullptr; const char* ref = nullptr; const char* alt = nullptr;
  std::size_t chr_n = 0, pos_n = 0, ref_n = 0, alt_n = 0;
};
inline bool variant_fields(const char* s, std::size_t len, bool bim, VariantFields& f) {
  const char* tok[6]; std::size_t tn[6];
  const int want = bim ? 6 : 5;
  const char* p = s; const char* const e = s + len;
  int k = 0;
  while (k < want) {
    while (p < e && (*p == '\t' || *p == ' ')) ++p;
    if (p == e) break;
    const char* b = p;
    while (p < e && *p != '\t' && *p != ' ') ++p;
    tok[k] = b; tn[k] = static_cast<std::size_t>(p - b); ++k;
  }
  if (k < want) return false;
  f.chr = tok[0]; f.chr_n = tn[0];
  if (bim) { f.pos = tok[3]; f.pos_n = tn[3]; f.alt = tok[4]; f.alt_n = tn[4]; f.ref = tok[5]; f.ref_n = tn[5]; }
  else { f.pos = tok[1]; f.pos_n = tn[1]; f.ref = tok[3]; f.ref_n = tn[3]; f.alt = tok[4]; f.alt_n = tn[4]; }
  return true;
}

inline char ascii_upper(char c) { return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : c; }

// The chromosome code of a label, the rule of R/chr.R: 1..22, 23 = X (X, 23, 25, PAR1, PAR2, XY), 24 = Y, 26 = MT
// (MT, M), with or without a leading "chr", any case; -1 for any other label.
inline int chr_code(const char* s, std::size_t n) {
  if (n >= 3 && ascii_upper(s[0]) == 'C' && ascii_upper(s[1]) == 'H' && ascii_upper(s[2]) == 'R') { s += 3; n -= 3; }
  if (n == 0 || n > 4) return -1;
  char u[4];
  for (std::size_t i = 0; i < n; ++i) u[i] = ascii_upper(s[i]);
  const auto is = [&](const char* w) { return std::strlen(w) == n && std::memcmp(u, w, n) == 0; };
  const auto digit = [](char c) { return c >= '0' && c <= '9'; };
  int code = -1;
  if (n <= 2 && digit(u[0]) && (n == 1 || digit(u[1]))) code = n == 1 ? u[0] - '0' : (u[0] - '0') * 10 + (u[1] - '0');
  else if (is("X") || is("PAR1") || is("PAR2") || is("XY")) code = 23;
  else if (is("Y")) code = 24;
  else if (is("MT") || is("M")) code = 26;
  if (code == 25) code = 23;
  return ((code >= 1 && code <= 24) || code == 26) ? code : -1;
}

// A base-pair position: digits only, below 2^31; -1 for anything else.
inline std::int32_t parse_bp(const char* s, std::size_t n) {
  if (n == 0 || n > 10) return -1;
  std::int64_t v = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (s[i] < '0' || s[i] > '9') return -1;
    v = v * 10 + (s[i] - '0');
  }
  return v > 2147483647LL ? -1 : static_cast<std::int32_t>(v);
}

// The number of alleles an ALT field names: the reference and one per comma-separated ALT allele.
inline int allele_count(const char* alt, std::size_t n) {
  int c = 2;
  for (std::size_t i = 0; i < n; ++i) c += alt[i] == ',';
  return c;
}

// The number of alleles of every variant line of a .pvar.
inline std::vector<int> allele_counts(const std::string& pvar) {
  VariantLines in(pvar, false);
  std::vector<int> ct;
  const char* s; std::size_t len;
  VariantFields f;
  while (in.next(s, len)) {
    if (len == 0) continue;
    if (!variant_fields(s, len, false, f))
      Rcpp::stop("Line %.0f of the variants of %s has fewer than 5 fields.", static_cast<double>(ct.size() + 1),
                 pvar.c_str());
    ct.push_back(allele_count(f.alt, f.alt_n));
  }
  return ct;
}

// The distinct allele strings of a table, in upper case, each with the code of its first appearance.
class AlleleCodes {
 public:
  AlleleCodes() { for (int& c : one_) c = -1; }
  std::int32_t code(const char* s, std::size_t n) {
    if (n == 1) {
      const unsigned char u = static_cast<unsigned char>(ascii_upper(s[0]));
      if (one_[u] < 0) { one_[u] = static_cast<std::int32_t>(text.size()); text.emplace_back(1, static_cast<char>(u)); }
      return one_[u];
    }
    key_.assign(s, n);
    for (char& c : key_) c = ascii_upper(c);
    const auto it = more_.find(key_);
    if (it != more_.end()) return it->second;
    const std::int32_t c = static_cast<std::int32_t>(text.size());
    text.push_back(key_);
    more_.emplace(key_, c);
    return c;
  }
  std::vector<std::string> text;

 private:
  std::int32_t one_[256];
  std::unordered_map<std::string, std::int32_t> more_;
  std::string key_;
};

// Every variant line of a .bim / .pvar, by file row: chromosome code and position (-1: none), the two alleles as
// codes into alleles (upper case), and for a .pvar the number of alleles (a .bim line has two: allele_ct is empty).
struct VariantTable {
  std::vector<std::int32_t> chr, bp, ref, alt;
  std::vector<int> allele_ct;
  std::vector<std::string> alleles;
  std::size_t size() const { return chr.size(); }
  bool biallelic(std::size_t r) const { return allele_ct.empty() || allele_ct[r] == 2; }
};
inline VariantTable read_variant_table(const std::string& path, bool bim) {
  VariantLines in(path, bim);
  VariantTable v;
  AlleleCodes codes;
  const char* s; std::size_t len;
  VariantFields f;
  while (in.next(s, len)) {
    if (len == 0) continue;
    if (!variant_fields(s, len, bim, f))
      Rcpp::stop("Line %.0f of the variants of %s has fewer than %d fields.", static_cast<double>(v.size() + 1),
                 path.c_str(), bim ? 6 : 5);
    v.chr.push_back(chr_code(f.chr, f.chr_n));
    v.bp.push_back(parse_bp(f.pos, f.pos_n));
    v.ref.push_back(codes.code(f.ref, f.ref_n));
    v.alt.push_back(codes.code(f.alt, f.alt_n));
    if (!bim) v.allele_ct.push_back(allele_count(f.alt, f.alt_n));
  }
  v.alleles = std::move(codes.text);
  return v;
}

// A number as the shortest text that reads back as the same double; NA for a missing value.
inline void put_number(std::string& o, double v) {
  if (std::isnan(v)) { o += "NA"; return; }
  if (std::isinf(v)) { o += v > 0 ? "Inf" : "-Inf"; return; }
  char b[32];
  o.append(b, std::to_chars(b, b + 32, v).ptr);
}
inline void put_count(std::string& o, long long v) {
  char b[24];
  o.append(b, std::to_chars(b, b + 24, v).ptr);
}

// The result of a scan of every variant of a file, as text: each line of the .bim / .pvar as it is, then the
// statistics of that variant, tab-separated; a variant without statistics has NA in their place. The first line names
// the columns: the "#CHROM" line of the .pvar, or CHR ID CM POS A1 A2 for a .bim (which has none), then the
// statistics. The blanks between the fields of a .bim line are written as tabs.
//
// The scan goes through the file in order. Before a group of threads computes, load() brings in the lines up to the
// group's last variant; every thread then writes the rows of its own variants into a text of its own (row()), and
// write() puts those texts into the file in order. finish() writes the lines after the last variant tested.
class ResultText {
 public:
  ResultText(const std::string& text_path, bool bim, const std::string& out_path,
             const std::vector<std::string>& statistics)
      : in_(text_path, bim), bim_(bim), n_stat_(static_cast<int>(statistics.size())), out_path_(out_path) {
    out_ = std::fopen(out_path.c_str(), "wb");
    if (!out_) Rcpp::stop("Cannot write " + out_path);
    std::setvbuf(out_, nullptr, _IOFBF, std::size_t(1) << 22);
    std::string h = bim ? std::string("CHR\tID\tCM\tPOS\tA1\tA2") : in_.header();
    for (const std::string& s : statistics) { h += '\t'; h += s; }
    h += '\n';
    put(h);
  }
  ~ResultText() {
    if (out_) { std::fclose(out_); std::remove(out_path_.c_str()); }  // not finished: no partial file is left
  }
  ResultText(const ResultText&) = delete;
  ResultText& operator=(const ResultText&) = delete;

  // Brings the lines of the file rows first() .. row_end - 1 into memory, in place of those held before.
  void load(long long row_end) {
    first_ = next_;
    text_.clear(); off_.clear();
    const char* s; std::size_t len;
    while (next_ < row_end) {
      if (!take(s, len))
        Rcpp::stop("%s has %.0f variant lines, fewer than the genotype file has variants.", in_.path().c_str(),
                   static_cast<double>(next_));
      off_.push_back(text_.size());
      append_line(text_, s, len);
      ++next_;
    }
    off_.push_back(text_.size());
  }
  long long first() const { return first_; }
  long long end() const { return next_; }
  // Appends the row of file row r (loaded) to o: its line, then n and the statistics v, or NA for each when v is null.
  void row(std::string& o, long long r, long long n, const double* v) const {
    const std::size_t i = static_cast<std::size_t>(r - first_);
    o.append(text_.data() + off_[i], off_[i + 1] - off_[i]);
    if (!v) {
      for (int k = 0; k < n_stat_; ++k) o += "\tNA";
    } else {
      o += '\t'; put_count(o, n);
      for (int k = 0; k + 1 < n_stat_; ++k) { o += '\t'; put_number(o, v[k]); }
    }
    o += '\n';
  }
  void write(const std::string& rows) { put(rows); }
  // The lines after the last variant tested (no statistics), the end of the file; returns the variant lines of the
  // .bim / .pvar, which must be m.
  long long finish(long long m) {
    std::string o;
    const char* s; std::size_t len;
    while (take(s, len)) {
      ++next_;
      append_line(o, s, len);
      for (int k = 0; k < n_stat_; ++k) o += "\tNA";
      o += '\n';
      if (o.size() > (std::size_t(1) << 20)) { put(o); o.clear(); }
    }
    put(o);
    if (next_ != m)
      Rcpp::stop("%s has %.0f variant lines, the genotype file %.0f variants.", in_.path().c_str(),
                 static_cast<double>(next_), static_cast<double>(m));
    const bool bad = std::fclose(out_) != 0;
    out_ = nullptr;
    if (bad) { std::remove(out_path_.c_str()); Rcpp::stop("Write error in " + out_path_); }
    return next_;
  }
  double bytes() const { return bytes_; }

 private:
  bool take(const char*& s, std::size_t& len) {  // the next variant line (empty lines are passed over)
    while (in_.next(s, len)) if (len) return true;
    return false;
  }
  void append_line(std::string& o, const char* s, std::size_t len) const {
    if (!bim_ || !std::memchr(s, ' ', len)) { o.append(s, len); return; }  // tab-separated already
    const char* p = s; const char* const e = s + len;
    bool any = false;
    while (p < e) {
      while (p < e && (*p == '\t' || *p == ' ')) ++p;
      if (p == e) break;
      const char* b = p;
      while (p < e && *p != '\t' && *p != ' ') ++p;
      if (any) o += '\t';
      o.append(b, static_cast<std::size_t>(p - b));
      any = true;
    }
  }
  void put(const std::string& s) {
    if (s.empty()) return;
    if (std::fwrite(s.data(), 1, s.size(), out_) != s.size()) Rcpp::stop("Write error in " + out_path_);
    bytes_ += static_cast<double>(s.size());
  }
  VariantLines in_;
  bool bim_;
  int n_stat_;
  std::string out_path_;
  std::FILE* out_ = nullptr;
  std::string text_;
  std::vector<std::size_t> off_;
  long long first_ = 0, next_ = 0;
  double bytes_ = 0;
};

}  // namespace plink4r

#endif
