// The variant text of a genotype file, read in compiled code: the .bim of a BED or the .pvar of a PGEN, line by line
// in file order. R holds none of it. Three layers:
//   TextLines      the lines of a text file, one after another
//   VariantLines   the variant lines of a .bim / .pvar, and what a field of such a line says (chr_code, parse_bp,
//                  variant_fields)
//   VariantTable   the variants of a file, read once (read_variant_table): a line that names K ALT alleles is K
//                  variants, each counting one of them, so every variant has one counted allele (A1) and one other (A2)
//   SampleTable    the samples of a .fam / .psam, read once (read_sample_table)
//   ResultText     the result of a scan as text: CHR, BP, A1, A2 of each variant scanned, then its statistics
// What a result says of a variant is what its line says: chromosome label and alleles are written as the file has
// them, never recoded. The chromosome code (chr_code) only matches a variant to a table that names chromosomes.
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

// The variant lines of a .bim (six columns, no header) or of a .pvar (its "##" lines and its "#CHROM" line are
// passed over).
class VariantLines {
 public:
  VariantLines(const std::string& path, bool bim) : in_(path) {
    if (bim) return;
    const char* s; std::size_t len;
    while (in_.next(s, len)) {
      if (len >= 2 && s[0] == '#' && s[1] == '#') continue;
      if (len >= 6 && std::memcmp(s, "#CHROM", 6) == 0) return;
      break;
    }
    Rcpp::stop("No #CHROM line before the variants of " + path);
  }
  const std::string& path() const { return in_.path(); }
  bool next(const char*& s, std::size_t& len) { return in_.next(s, len); }

 private:
  TextLines in_;
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

// The distinct texts of a field, as the file writes them, each with the code of its first appearance.
class TextCodes {
 public:
  TextCodes() { for (int& c : one_) c = -1; }
  std::int32_t code(const char* s, std::size_t n) {
    if (n == 1) {
      const unsigned char u = static_cast<unsigned char>(s[0]);
      if (one_[u] < 0) { one_[u] = static_cast<std::int32_t>(text.size()); text.emplace_back(1, s[0]); }
      return one_[u];
    }
    key_.assign(s, n);
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

// The variants of a .bim / .pvar, in file order. A .pvar line that names K ALT alleles (comma-separated) is K
// variants, one per ALT allele, in the order the line names them: variant v counts the alt[v]-th ALT of file line
// line[v] (A1) against the reference (A2); a sample carrying another ALT of the line counts 0 of this one. A .bim
// line is one variant: A1 = column 5, the counted allele, A2 = column 6. chr, a1, a2: codes into the distinct texts
// of the file (chromosomes, alleles: as the file writes them); bp: the position (-1: not a position);
// chromosome_code: the chromosome code of each label (chr_code; -1: none). line_alleles: the number of alleles of
// every file line (a PGEN needs it before any genotype is read).
struct VariantTable {
  std::vector<std::int32_t> chr, bp, a1, a2, line, alt;
  std::vector<int> line_alleles;
  std::vector<std::string> chromosomes, alleles;
  std::vector<std::int32_t> chromosome_code;
  std::size_t size() const { return chr.size(); }
  std::size_t lines() const { return line_alleles.size(); }
  std::int32_t code(std::size_t v) const { return chromosome_code[chr[v]]; }
};
inline VariantTable read_variant_table(const std::string& path, bool bim) {
  VariantLines in(path, bim);
  VariantTable v;
  TextCodes chromosomes, alleles;
  std::int32_t c = -1;  // the label of the line before: the lines of a chromosome follow one another
  const char* s; std::size_t len;
  VariantFields f;
  while (in.next(s, len)) {
    if (len == 0) continue;
    const std::int32_t line = static_cast<std::int32_t>(v.lines());
    if (!variant_fields(s, len, bim, f))
      Rcpp::stop("Line %.0f of the variants of %s has fewer than %d fields.", static_cast<double>(line + 1),
                 path.c_str(), bim ? 6 : 5);
    if (c < 0 || chromosomes.text[c].size() != f.chr_n || std::memcmp(chromosomes.text[c].data(), f.chr, f.chr_n) != 0)
      c = chromosomes.code(f.chr, f.chr_n);
    const std::int32_t bp = parse_bp(f.pos, f.pos_n), ref = alleles.code(f.ref, f.ref_n);
    int k = 0;
    for (const char* b = f.alt; ; ) {  // the ALT alleles of the line: one for a .bim line
      const char* e = bim ? nullptr : static_cast<const char*>(std::memchr(b, ',', f.alt + f.alt_n - b));
      if (!e) e = f.alt + f.alt_n;
      v.chr.push_back(c); v.bp.push_back(bp); v.a2.push_back(ref);
      v.a1.push_back(alleles.code(b, static_cast<std::size_t>(e - b)));
      v.line.push_back(line); v.alt.push_back(++k);
      if (e == f.alt + f.alt_n) break;
      b = e + 1;
    }
    v.line_alleles.push_back(k + 1);
  }
  v.chromosomes = std::move(chromosomes.text);
  v.alleles = std::move(alleles.text);
  for (const std::string& label : v.chromosomes) v.chromosome_code.push_back(chr_code(label.data(), label.size()));
  return v;
}

// The samples of a .fam (columns 1 and 2: family and individual ID, no header) or of a .psam (the columns named
// IID or #IID, and FID or #FID in its header line; fid "0" without that column), in file order.
struct SampleTable {
  std::vector<std::string> fid, iid;
  std::size_t size() const { return iid.size(); }
};
inline SampleTable read_sample_table(const std::string& path, bool fam) {
  TextLines in(path);
  SampleTable t;
  const char* s; std::size_t len;
  std::vector<std::pair<const char*, std::size_t>> tok;
  const auto split = [&](bool tabs_only) {
    tok.clear();
    const char* p = s; const char* const e = s + len;
    if (tabs_only) {  // a .psam: tab-separated, a field may be empty
      for (;;) {
        const char* q = static_cast<const char*>(std::memchr(p, '\t', e - p));
        if (!q) q = e;
        tok.emplace_back(p, static_cast<std::size_t>(q - p));
        if (q == e) break;
        p = q + 1;
      }
      return;
    }
    while (p < e) {
      while (p < e && (*p == '\t' || *p == ' ')) ++p;
      if (p == e) break;
      const char* b = p;
      while (p < e && *p != '\t' && *p != ' ') ++p;
      tok.emplace_back(b, static_cast<std::size_t>(p - b));
    }
  };
  int fid_col = fam ? 0 : -1, iid_col = fam ? 1 : -1;
  bool header = !fam, tabs = false;
  while (in.next(s, len)) {
    if (len == 0) continue;
    if (header) {
      if (len >= 2 && s[0] == '#' && s[1] == '#') continue;
      tabs = std::memchr(s, '\t', len) != nullptr;
      split(tabs);
      for (std::size_t c = 0; c < tok.size(); ++c) {
        const std::string name(tok[c].first, tok[c].second);
        if (name == "#IID" || name == "IID") { if (iid_col < 0) iid_col = static_cast<int>(c); }
        else if (name == "#FID" || name == "FID") { if (fid_col < 0) fid_col = static_cast<int>(c); }
      }
      if (iid_col < 0) Rcpp::stop(path + " has no IID column");
      header = false;
      continue;
    }
    split(tabs);
    if (static_cast<int>(tok.size()) <= std::max(iid_col, fid_col))
      Rcpp::stop("Line %.0f of the samples of %s has too few fields.", static_cast<double>(t.size() + 1), path.c_str());
    t.iid.emplace_back(tok[iid_col].first, tok[iid_col].second);
    if (fid_col >= 0) t.fid.emplace_back(tok[fid_col].first, tok[fid_col].second); else t.fid.emplace_back("0");
  }
  return t;
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

// The result of a scan as text, tab-separated: one row per variant scanned, in the order of the scan. A row is CHR,
// BP, A1, A2 of the variant (chromosome label and alleles as its line of the .bim / .pvar writes them, nothing
// recoded; A1 is the counted allele; NA for a position that is not one), then its statistics; a variant without
// statistics has NA in their place. The first line names the columns. Every thread writes the rows of its own
// variants into a text of its own (row()), and write() puts those texts into the file in order; finish() closes it.
class ResultText {
 public:
  // statistics: the names of the columns after A2; the second of them is a count (row()).
  ResultText(const VariantTable& variants, const std::string& out_path, const std::vector<std::string>& statistics)
      : vt_(variants), n_stat_(static_cast<int>(statistics.size())), out_path_(out_path) {
    out_ = std::fopen(out_path.c_str(), "wb");
    if (!out_) Rcpp::stop("Cannot write " + out_path);
    std::setvbuf(out_, nullptr, _IOFBF, std::size_t(1) << 22);
    std::string h = "CHR\tBP\tA1\tA2";
    for (const std::string& s : statistics) { h += '\t'; h += s; }
    h += '\n';
    put(h);
  }
  ~ResultText() {
    if (out_) { std::fclose(out_); std::remove(out_path_.c_str()); }  // not finished: no partial file is left
  }
  ResultText(const ResultText&) = delete;
  ResultText& operator=(const ResultText&) = delete;

  // Appends the row of variant r of the table to o: its four fields, then the statistics: v[0], the count n, v[1],
  // v[2], ...
  void row(std::string& o, std::size_t r, long long n, const double* v) const {
    o += vt_.chromosomes[vt_.chr[r]]; o += '\t';
    if (vt_.bp[r] < 0) o += "NA"; else put_count(o, vt_.bp[r]);
    o += '\t'; o += vt_.alleles[vt_.a1[r]];
    o += '\t'; o += vt_.alleles[vt_.a2[r]];
    o += '\t'; put_number(o, v[0]);
    o += '\t'; put_count(o, n);
    for (int k = 1; k + 1 < n_stat_; ++k) { o += '\t'; put_number(o, v[k]); }
    o += '\n';
  }
  void write(const std::string& rows) { put(rows); }
  void finish() {
    const bool bad = std::fclose(out_) != 0;
    out_ = nullptr;
    if (bad) { std::remove(out_path_.c_str()); Rcpp::stop("Write error in " + out_path_); }
  }
  double bytes() const { return bytes_; }

 private:
  void put(const std::string& s) {
    if (s.empty()) return;
    if (std::fwrite(s.data(), 1, s.size(), out_) != s.size()) Rcpp::stop("Write error in " + out_path_);
    bytes_ += static_cast<double>(s.size());
  }
  const VariantTable& vt_;
  int n_stat_;
  std::string out_path_;
  std::FILE* out_ = nullptr;
  double bytes_ = 0;
};

}  // namespace plink4r

#endif
