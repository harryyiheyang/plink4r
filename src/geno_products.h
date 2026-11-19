// The product kernels of a geno object (geno_products.cpp) for the C-callable entry points (geno_api.cpp).
// All write into caller buffers, column-major; cols = nullptr means every column.
#ifndef PLINK4R_GENO_PRODUCTS_H
#define PLINK4R_GENO_PRODUCTS_H

#include "geno_store.h"

namespace plink4r {
namespace gstore {

void xtx_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads);
void xtwx_float_into(const Geno& g, const double* w, const int* cols, std::size_t pc, double* o, int threads);
void xtwx_bin_into(const Geno& g, const double* w, const int* cols, std::size_t pc, double* o, double* w_bin,
                   int threads);
void xtwx_std(const Geno& g, const double* w, const int* cols, std::size_t pc, const double* center,
              const double* div, double* o, int threads);
void xtm_into(const Geno& g, const double* m, std::size_t k, double* o, int threads);
void mx_into(const Geno& g, const double* b, std::size_t k, double* o, int threads);
void dense_into(const Geno& g, const int* rows, std::size_t nr, const int* cols, std::size_t nc, double* o,
                int threads);
void std_into(const double* k, std::size_t r, std::size_t c, const double* u, const double* v, double s,
              const double* da, const double* db, double* o, int threads);

}  // namespace gstore
}  // namespace plink4r

#endif
