// Matching of variant tables by the quadruple (chromosome code, position, A1, A2), free of Rcpp types.
//
// A query row matches a reference row with the same chromosome code (R/chr.R) and position whose allele pair is
// the same (sign +1: the query counts the reference's A1) or swapped (sign -1). Alleles are compared as given: pass
// them through upper_all() first to compare without case. No strand flipping. Each reference row is used once and
// each query row once: query rows are taken in order, the same orientation before the swapped one. A query row whose
// candidates are all used already is a duplicate (row -2); a row without any candidate is unmatched (row -1). A
// position equal to `na` (R's NA_integer_) or a chromosome code <= 0 never matches. quad_dup() flags the rows of one
// table whose variant (position and unordered allele pair) occurs more than once, every copy.
#ifndef PLINK4R_QUAD_MATCH_H
#define PLINK4R_QUAD_MATCH_H

#include <cctype>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace plink4r {

struct Quad {
  const std::vector<int>* chr;
  const std::vector<int>* bp;
  const std::vector<std::string>* a1;
  const std::vector<std::string>* a2;
  std::size_t size() const { return chr->size(); }
};

struct QuadHit {
  std::vector<int> row;   // reference row (0-based), -1 unmatched, -2 duplicate
  std::vector<int> sign;  // +1 same pair, -1 swapped, 0 otherwise
  int n = 0;              // matched rows
};

inline std::string upper(const std::string& s) {
  std::string u(s);
  for (char& c : u) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return u;
}

inline std::vector<std::string> upper_all(const std::vector<std::string>& x) {
  std::vector<std::string> u(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) u[i] = upper(x[i]);
  return u;
}

inline std::vector<char> quad_dup(const Quad& q, int na) {
  std::unordered_map<std::string, int> first;
  first.reserve(q.size());
  std::vector<char> dup(q.size(), 0);
  for (std::size_t j = 0; j < q.size(); ++j) {
    const int c = (*q.chr)[j], p = (*q.bp)[j];
    if (c <= 0 || p == na) continue;
    const std::string& a = (*q.a1)[j];
    const std::string& b = (*q.a2)[j];
    const std::string key = std::to_string(c) + ':' + std::to_string(p) + ':' + (a < b ? a + '/' + b : b + '/' + a);
    const auto it = first.emplace(key, static_cast<int>(j));
    if (!it.second) { dup[j] = 1; dup[it.first->second] = 1; }
  }
  return dup;
}

inline QuadHit quad_match(const Quad& ref, const Quad& qry, int na) {
  auto key = [](int c, int p) { return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(c)) << 32) | static_cast<std::uint32_t>(p); };
  std::unordered_map<std::uint64_t, std::vector<int>> pos;
  pos.reserve(ref.size());
  for (std::size_t i = 0; i < ref.size(); ++i) {
    const int c = (*ref.chr)[i], p = (*ref.bp)[i];
    if (c > 0 && p != na) pos[key(c, p)].push_back(static_cast<int>(i));
  }
  QuadHit h;
  h.row.assign(qry.size(), -1);
  h.sign.assign(qry.size(), 0);
  std::vector<char> used(ref.size(), 0);
  for (std::size_t j = 0; j < qry.size(); ++j) {
    const int c = (*qry.chr)[j], p = (*qry.bp)[j];
    if (c <= 0 || p == na) continue;
    const auto it = pos.find(key(c, p));
    if (it == pos.end()) continue;
    const std::string& b1 = (*qry.a1)[j];
    const std::string& b2 = (*qry.a2)[j];
    bool cand = false;
    for (int pass = 0; pass < 2 && h.row[j] < 0; ++pass) {
      for (int i : it->second) {
        const std::string& r1 = (*ref.a1)[i];
        const std::string& r2 = (*ref.a2)[i];
        const bool ok = pass == 0 ? (r1 == b1 && r2 == b2) : (r1 == b2 && r2 == b1);
        if (!ok) continue;
        cand = true;
        if (used[i]) continue;
        used[i] = 1;
        h.row[j] = i;
        h.sign[j] = pass == 0 ? 1 : -1;
        ++h.n;
        break;
      }
    }
    if (h.row[j] < 0 && cand) h.row[j] = -2;
  }
  return h;
}

// The reference has two position columns (builds); the query is matched on each and the build that matches more
// rows is kept (the first on a tie). build_out: 0 or 1.
inline QuadHit quad_match_best(const Quad& ref0, const Quad& ref1, const Quad& qry, int na, int& build_out) {
  QuadHit h0 = quad_match(ref0, qry, na);
  QuadHit h1 = quad_match(ref1, qry, na);
  build_out = h1.n > h0.n ? 1 : 0;
  return build_out ? h1 : h0;
}

}  // namespace plink4r
#endif
