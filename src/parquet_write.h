// A Parquet file written row group by row group: the package's Parquet writer for results that go to disk as they are
// computed, with no table held in memory. A flat schema of required columns (no nulls, no nesting); pages are not
// compressed; one data page per column and row group. Encodings, chosen by the caller per column:
//   float_plain    FLOAT, PLAIN (4 bytes per value)
//   int32_delta    INT32, DELTA_BINARY_PACKED: for values that change slowly, e.g. positions in order
//   int32_dict     INT32, RLE_DICTIONARY: for few distinct values in a row group (one value: a few bytes for the group)
//   string_dict    strings given as codes into a table of levels, RLE_DICTIONARY (the dictionary is per row group)
//   string_const   one string for the whole row group
// INT32 and string column chunks carry their minimum and maximum, so a reader skips the row groups a filter excludes.
// A file that is not closed (an error or an interrupt on the way) is removed.
// Format: parquet.apache.org/docs/file-format and parquet.thrift (footer and page headers: Thrift compact protocol).
#ifndef PLINK4R_PARQUET_WRITE_H
#define PLINK4R_PARQUET_WRITE_H

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace plink4r {
namespace pq {

using Bytes = std::vector<unsigned char>;

inline void put_varint(Bytes& b, std::uint64_t v) {
  while (v >= 0x80) { b.push_back(static_cast<unsigned char>((v & 0x7f) | 0x80)); v >>= 7; }
  b.push_back(static_cast<unsigned char>(v));
}
inline std::uint64_t zigzag(std::int64_t v) {
  return (static_cast<std::uint64_t>(v) << 1) ^ static_cast<std::uint64_t>(v >> 63);
}
// Bits needed for v (0 for 0).
inline int bits_of(std::uint32_t v) {
  int n = 0;
  while (v) { ++n; v >>= 1; }
  return n;
}

// Values of a fixed number of bits each, least significant bit first.
struct BitPacker {
  Bytes& out;
  std::uint64_t acc = 0;
  int nbits = 0;
  explicit BitPacker(Bytes& o) : out(o) {}
  void put(std::uint32_t v, int width) {
    acc |= static_cast<std::uint64_t>(v) << nbits;
    nbits += width;
    while (nbits >= 8) { out.push_back(static_cast<unsigned char>(acc & 0xff)); acc >>= 8; nbits -= 8; }
  }
  void flush() {
    if (nbits > 0) { out.push_back(static_cast<unsigned char>(acc & 0xff)); acc = 0; nbits = 0; }
  }
};

// Thrift compact protocol: a struct is open on construction; fields are written in increasing id.
struct Thrift {
  enum { kI32 = 5, kI64 = 6, kBinary = 8, kList = 9, kStruct = 12 };
  Bytes b;
  std::vector<int> last{0};  // last field id of every open struct
  void field(int id, int type) {
    const int delta = id - last.back();
    if (delta > 0 && delta <= 15) b.push_back(static_cast<unsigned char>((delta << 4) | type));
    else { b.push_back(static_cast<unsigned char>(type)); put_varint(b, zigzag(id)); }
    last.back() = id;
  }
  void i32(int id, std::int64_t v) { field(id, kI32); put_varint(b, zigzag(v)); }
  void i64(int id, std::int64_t v) { field(id, kI64); put_varint(b, zigzag(v)); }
  void raw(const void* p, std::size_t n) {
    put_varint(b, n);
    const unsigned char* c = static_cast<const unsigned char*>(p);
    b.insert(b.end(), c, c + n);
  }
  void binary(int id, const void* p, std::size_t n) { field(id, kBinary); raw(p, n); }
  void str(int id, const std::string& s) { binary(id, s.data(), s.size()); }
  void list(int id, int elem_type, std::size_t n) {
    field(id, kList);
    if (n < 15) b.push_back(static_cast<unsigned char>((n << 4) | elem_type));
    else { b.push_back(static_cast<unsigned char>(0xF0 | elem_type)); put_varint(b, n); }
  }
  void elem_i32(std::int64_t v) { put_varint(b, zigzag(v)); }
  void elem_str(const std::string& s) { raw(s.data(), s.size()); }
  void begin(int id) { field(id, kStruct); last.push_back(0); }  // a struct field
  void begin_elem() { last.push_back(0); }                        // a struct of a list
  void end() { b.push_back(0); last.pop_back(); }
};

// Dictionary indices of a data page: the bit width, then runs of the RLE / bit-packed hybrid. A run of 8 or more
// equal indices is one RLE run; the rest go bit-packed in groups of 8.
inline void put_indices(Bytes& out, const std::uint32_t* idx, std::size_t n, std::uint32_t n_dict) {
  const int width = std::max(1, bits_of(n_dict > 0 ? n_dict - 1 : 0));
  const int vbytes = (width + 7) / 8;
  out.push_back(static_cast<unsigned char>(width));
  auto run_at = [&](std::size_t i) {  // length of the run of equal indices that starts at i
    std::size_t r = 1;
    while (i + r < n && idx[i + r] == idx[i]) ++r;
    return r;
  };
  std::size_t i = 0;
  while (i < n) {
    const std::size_t r = run_at(i);
    if (r >= 8) {
      put_varint(out, static_cast<std::uint64_t>(r) << 1);
      for (int k = 0; k < vbytes; ++k) out.push_back(static_cast<unsigned char>((idx[i] >> (8 * k)) & 0xff));
      i += r;
      continue;
    }
    const std::size_t start = i;
    for (;;) {  // whole groups of 8, until the end or a run of 8 equal indices starts at a group boundary
      i = std::min(n, i + 8);
      if (i >= n || run_at(i) >= 8) break;
    }
    const std::size_t groups = (i - start + 7) / 8;
    put_varint(out, (static_cast<std::uint64_t>(groups) << 1) | 1);
    BitPacker bp(out);
    for (std::size_t k = 0; k < groups * 8; ++k) bp.put(start + k < n ? idx[start + k] : 0, width);
    bp.flush();
  }
}

// DELTA_BINARY_PACKED of INT32: blocks of 128 deltas in 4 miniblocks of 32, each miniblock at its own bit width.
inline void put_delta_i32(Bytes& out, const std::int32_t* v, std::size_t n) {
  constexpr std::size_t kBlock = 128, kMini = 4, kPer = kBlock / kMini;
  put_varint(out, kBlock); put_varint(out, kMini); put_varint(out, n);
  put_varint(out, zigzag(n ? v[0] : 0));
  std::int32_t delta[kBlock];
  std::uint32_t rel[kBlock];
  for (std::size_t i = 1; i < n;) {
    const std::size_t m = std::min(kBlock, n - i);
    std::int32_t low = INT32_MAX;
    for (std::size_t c = 0; c < m; ++c) {
      delta[c] = static_cast<std::int32_t>(static_cast<std::uint32_t>(v[i + c]) - static_cast<std::uint32_t>(v[i + c - 1]));
      low = std::min(low, delta[c]);
    }
    for (std::size_t c = 0; c < m; ++c) rel[c] = static_cast<std::uint32_t>(delta[c]) - static_cast<std::uint32_t>(low);
    put_varint(out, zigzag(low));
    int width[kMini];
    for (std::size_t mb = 0; mb < kMini; ++mb) {
      std::uint32_t top = 0;
      for (std::size_t c = mb * kPer; c < std::min(m, (mb + 1) * kPer); ++c) top = std::max(top, rel[c]);
      width[mb] = bits_of(top);
      out.push_back(static_cast<unsigned char>(width[mb]));
    }
    for (std::size_t mb = 0; mb * kPer < m; ++mb) {  // the last miniblock is filled with zeros
      BitPacker bp(out);
      for (std::size_t c = mb * kPer; c < (mb + 1) * kPer; ++c) bp.put(c < m ? rel[c] : 0, width[mb]);
      bp.flush();
    }
    i += m;
  }
}

}  // namespace pq

class ParquetWriter {
 public:
  enum Type { kInt32 = 1, kFloat = 4, kString = 6 };  // parquet physical types (kString: BYTE_ARRAY, UTF8)
  struct Column { std::string name; Type type; };

  // meta: key / value pairs of the file. created_by: "<application> version <x.y.z>".
  ParquetWriter(const std::string& path, std::vector<Column> columns,
                std::vector<std::pair<std::string, std::string>> meta, std::string created_by)
      : path_(path), cols_(std::move(columns)), meta_(std::move(meta)), created_by_(std::move(created_by)),
        out_(path, std::ios::binary | std::ios::trunc) {
    if (!out_) throw std::runtime_error("Cannot write " + path);
    out_.write("PAR1", 4);
    pos_ = 4;
  }
  ~ParquetWriter() {
    if (!closed_) { out_.close(); std::remove(path_.c_str()); }
  }
  ParquetWriter(const ParquetWriter&) = delete;
  ParquetWriter& operator=(const ParquetWriter&) = delete;

  // A row group of `rows` rows: then one call per column, in the order of the schema, each with `rows` values.
  void begin_group(std::size_t rows) {
    if (rows == 0 || rows > static_cast<std::size_t>(kMaxRows)) throw std::runtime_error("Parquet row group size");
    groups_.push_back(Group());
    groups_.back().rows = rows;
    rows_ = rows;
  }
  void float_plain(const float* v) {
    Chunk c = start(kFloat);
    c.encodings = {kPlain};
    page(c, false, kPlain, rows_, v, rows_ * sizeof(float));
    groups_.back().chunks.push_back(c);
  }
  void int32_delta(const std::int32_t* v) {
    Chunk c = start(kInt32);
    c.encodings = {kDelta};
    buf_.clear();
    pq::put_delta_i32(buf_, v, rows_);
    page(c, false, kDelta, rows_, buf_.data(), buf_.size());
    const auto mm = std::minmax_element(v, v + rows_);
    stats(c, le32(*mm.first), le32(*mm.second));
    groups_.back().chunks.push_back(c);
  }
  void int32_dict(const std::int32_t* v) {
    Chunk c = start(kInt32);
    uniq_.assign(v, v + rows_);
    std::sort(uniq_.begin(), uniq_.end());
    uniq_.erase(std::unique(uniq_.begin(), uniq_.end()), uniq_.end());
    idx_.resize(rows_);
    for (std::size_t i = 0; i < rows_; ++i)
      idx_[i] = static_cast<std::uint32_t>(std::lower_bound(uniq_.begin(), uniq_.end(), v[i]) - uniq_.begin());
    buf_.clear();
    for (std::int32_t u : uniq_) { const std::string s = le32(u); buf_.insert(buf_.end(), s.begin(), s.end()); }
    dictionary(c, uniq_.size());
    stats(c, le32(uniq_.front()), le32(uniq_.back()));
    groups_.back().chunks.push_back(c);
  }
  void string_dict(const std::int32_t* code, const std::vector<std::string>& levels) {
    Chunk c = start(kString);
    uniq_.assign(code, code + rows_);
    std::sort(uniq_.begin(), uniq_.end());
    uniq_.erase(std::unique(uniq_.begin(), uniq_.end()), uniq_.end());
    if (uniq_.front() < 0 || static_cast<std::size_t>(uniq_.back()) >= levels.size())
      throw std::runtime_error("Parquet string code out of range");
    idx_.resize(rows_);
    for (std::size_t i = 0; i < rows_; ++i)
      idx_[i] = static_cast<std::uint32_t>(std::lower_bound(uniq_.begin(), uniq_.end(), code[i]) - uniq_.begin());
    buf_.clear();
    const std::string* lo = &levels[uniq_.front()];
    const std::string* hi = lo;
    for (std::int32_t u : uniq_) {
      const std::string& s = levels[u];
      put_string(buf_, s);
      if (s < *lo) lo = &s;
      if (*hi < s) hi = &s;
    }
    dictionary(c, uniq_.size());
    stats(c, *lo, *hi);
    groups_.back().chunks.push_back(c);
  }
  void string_const(const std::string& s) {
    Chunk c = start(kString);
    idx_.assign(rows_, 0u);
    buf_.clear();
    put_string(buf_, s);
    dictionary(c, 1);
    stats(c, s, s);
    groups_.back().chunks.push_back(c);
  }
  void end_group() {
    Group& g = groups_.back();
    if (g.chunks.size() != cols_.size()) throw std::runtime_error("Parquet row group lacks columns");
    for (const Chunk& c : g.chunks) g.bytes += c.size;
    total_rows_ += static_cast<std::int64_t>(g.rows);
    if (!out_) throw std::runtime_error("Write failed: " + path_);
  }
  // Writes the footer; the file is complete after this.
  void close() {
    using pq::Thrift;
    Thrift t;
    t.i32(1, 2);
    t.list(2, Thrift::kStruct, cols_.size() + 1);
    t.begin_elem(); t.str(4, "schema"); t.i32(5, static_cast<std::int64_t>(cols_.size())); t.end();
    for (const Column& c : cols_) {
      t.begin_elem();
      t.i32(1, c.type); t.i32(3, 0); t.str(4, c.name);
      if (c.type == kString) { t.i32(6, 0); t.begin(10); t.begin(1); t.end(); t.end(); }  // UTF8; logical type STRING
      t.end();
    }
    t.i64(3, total_rows_);
    t.list(4, Thrift::kStruct, groups_.size());
    for (const Group& g : groups_) {
      t.begin_elem();
      t.list(1, Thrift::kStruct, g.chunks.size());
      for (std::size_t i = 0; i < g.chunks.size(); ++i) {
        const Chunk& c = g.chunks[i];
        t.begin_elem();
        t.i64(2, c.first);
        t.begin(3);
        t.i32(1, c.type);
        t.list(2, Thrift::kI32, c.encodings.size());
        for (int e : c.encodings) t.elem_i32(e);
        t.list(3, Thrift::kBinary, 1); t.elem_str(cols_[i].name);
        t.i32(4, 0);  // not compressed
        t.i64(5, static_cast<std::int64_t>(g.rows));
        t.i64(6, c.size); t.i64(7, c.size);
        t.i64(9, c.data);
        if (c.dict >= 0) t.i64(11, c.dict);
        if (c.has_stats) {
          t.begin(12); t.i64(3, 0); t.binary(5, c.max.data(), c.max.size()); t.binary(6, c.min.data(), c.min.size());
          t.end();
        }
        t.end();
        t.end();
      }
      t.i64(2, g.bytes); t.i64(3, static_cast<std::int64_t>(g.rows));
      t.i64(5, g.chunks.front().first); t.i64(6, g.bytes);
      t.end();
    }
    t.list(5, Thrift::kStruct, meta_.size());
    for (const auto& kv : meta_) { t.begin_elem(); t.str(1, kv.first); t.str(2, kv.second); t.end(); }
    t.str(6, created_by_);
    t.list(7, Thrift::kStruct, cols_.size());
    for (std::size_t i = 0; i < cols_.size(); ++i) { t.begin_elem(); t.begin(1); t.end(); t.end(); }  // type order
    t.end();
    if (t.b.size() > static_cast<std::size_t>(INT32_MAX)) throw std::runtime_error("Parquet footer too large");
    write(t.b.data(), t.b.size());
    const std::string len = le32(static_cast<std::int32_t>(t.b.size()));
    write(len.data(), 4);
    write("PAR1", 4);
    out_.close();
    if (out_.fail()) throw std::runtime_error("Write failed: " + path_);
    closed_ = true;
  }
  std::uint64_t bytes() const { return pos_; }
  std::int64_t rows() const { return total_rows_; }
  std::size_t groups() const { return groups_.size(); }

 private:
  enum { kPlain = 0, kRle = 3, kDelta = 5, kRleDictionary = 8 };
  static constexpr std::int64_t kMaxRows = 400000000;  // keeps every page below the 2^31 - 1 bytes of its header
  struct Chunk {
    int type = 0;
    std::vector<int> encodings;
    std::int64_t size = 0, first = 0, data = 0, dict = -1;
    bool has_stats = false;
    std::string min, max;
  };
  struct Group {
    std::vector<Chunk> chunks;
    std::size_t rows = 0;
    std::int64_t bytes = 0;
  };

  static std::string le32(std::int32_t v) {
    std::string s(4, '\0');
    const std::uint32_t u = static_cast<std::uint32_t>(v);
    for (int k = 0; k < 4; ++k) s[k] = static_cast<char>((u >> (8 * k)) & 0xff);
    return s;
  }
  static void put_string(pq::Bytes& b, const std::string& s) {
    const std::string len = le32(static_cast<std::int32_t>(s.size()));
    b.insert(b.end(), len.begin(), len.end());
    b.insert(b.end(), s.begin(), s.end());
  }
  void write(const void* p, std::size_t n) {
    out_.write(static_cast<const char*>(p), static_cast<std::streamsize>(n));
    pos_ += n;
  }
  Chunk start(Type type) {
    if (groups_.empty()) throw std::runtime_error("Parquet column outside a row group");
    const std::size_t i = groups_.back().chunks.size();
    if (i >= cols_.size() || cols_[i].type != type) throw std::runtime_error("Parquet column out of order");
    Chunk c;
    c.type = type;
    c.first = static_cast<std::int64_t>(pos_);
    return c;
  }
  void stats(Chunk& c, const std::string& lo, const std::string& hi) { c.has_stats = true; c.min = lo; c.max = hi; }
  // A page: its header, then `size` bytes. dict: a dictionary page (PLAIN), else a data page of the encoding.
  void page(Chunk& c, bool dict, int encoding, std::size_t n_values, const void* bytes, std::size_t size) {
    if (size > static_cast<std::size_t>(INT32_MAX)) throw std::runtime_error("Parquet page too large");
    pq::Thrift t;
    t.i32(1, dict ? 2 : 0);
    t.i32(2, static_cast<std::int64_t>(size)); t.i32(3, static_cast<std::int64_t>(size));
    if (dict) { t.begin(7); t.i32(1, static_cast<std::int64_t>(n_values)); t.i32(2, kPlain); t.end(); }
    else {
      t.begin(5); t.i32(1, static_cast<std::int64_t>(n_values)); t.i32(2, encoding); t.i32(3, kRle); t.i32(4, kRle);
      t.end();
    }
    t.end();
    if (dict) c.dict = static_cast<std::int64_t>(pos_); else c.data = static_cast<std::int64_t>(pos_);
    write(t.b.data(), t.b.size());
    write(bytes, size);
    c.size += static_cast<std::int64_t>(t.b.size() + size);
  }
  // The dictionary page held in buf_ (n_dict values), then the data page of idx_.
  void dictionary(Chunk& c, std::size_t n_dict) {
    c.encodings = {kPlain, kRleDictionary};
    page(c, true, kPlain, n_dict, buf_.data(), buf_.size());
    buf_.clear();
    pq::put_indices(buf_, idx_.data(), rows_, static_cast<std::uint32_t>(n_dict));
    page(c, false, kRleDictionary, rows_, buf_.data(), buf_.size());
  }

  std::string path_;
  std::vector<Column> cols_;
  std::vector<std::pair<std::string, std::string>> meta_;
  std::string created_by_;
  std::ofstream out_;
  std::uint64_t pos_ = 0;
  std::int64_t total_rows_ = 0;
  std::size_t rows_ = 0;
  bool closed_ = false;
  std::vector<Group> groups_;
  pq::Bytes buf_;
  std::vector<std::int32_t> uniq_;
  std::vector<std::uint32_t> idx_;
};

}  // namespace plink4r

#endif
