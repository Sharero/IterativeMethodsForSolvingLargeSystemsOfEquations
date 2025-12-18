#include "../include/symmetric_matrix.h"

#include <omp.h>

#include <algorithm>
#include <cstddef>

static constexpr bool IGNORE_OMP = false;
static constexpr size_t NUM_THREADS = 8;
static constexpr size_t MIN_SIZE = 1000;

#define VEC_SCHEDULE schedule(static)

std::optional<std::complex<double>> SymmetricMatrix::tryGetEntry(
    size_t row, size_t col) const {
    size_t size = this->size();

    if (row == col) {
        int blkBeg = idi[row];
        int blkSize = idi[row + 1] - blkBeg;

        double real = di[blkBeg];
        double imag = (blkSize == 2) ? di[blkBeg + 1] : 0.0;

        return std::complex<double>{real, imag};
    }

    if (col > row) {
        std::swap(row, col);
    }

    size_t start = static_cast<size_t>(ig[row]);
    size_t end = static_cast<size_t>(ig[row + 1]);

    if (start >= end) return std::nullopt;

    int target = static_cast<int>(col);
    auto it_begin = jg.begin() + static_cast<std::ptrdiff_t>(start);
    auto it_end = jg.begin() + static_cast<std::ptrdiff_t>(end);

    auto it = std::lower_bound(it_begin, it_end, target);

    if (it == it_end || *it != target) {
        return std::nullopt;
    }

    size_t el_id = static_cast<size_t>(std::distance(jg.begin(), it));

    int blkBeg = igg[el_id];
    int blkSize = igg[el_id + 1] - blkBeg;

    double real = gg[blkBeg];
    double imag = (blkSize == 2) ? gg[blkBeg + 1] : 0.0;

    return std::complex<double>{real, imag};
}

void SymmetricMatrix::multiplyMatrixBlockToVector(const CSpan& matSpan,
                                                  const CSpan& vecSpan,
                                                  Span resSpan) {
    auto mat = matSpan.begin();
    auto vec = vecSpan.begin();
    auto res = resSpan.begin();

    auto y0 = mat[0] * vec[0];
    auto y1 = mat[0] * vec[1];

    if (matSpan.size == 2) {
        y0 -= mat[1] * vec[1];
        y1 += mat[1] * vec[0];
    }

    res[0] += y0;
    res[1] += y1;
}

void SymmetricMatrix::multiplyMatrixByVector(const SymmetricMatrix& mat,
                                             const Vector& x, Vector* res) {
    res->clear();

    size_t size = mat.size();
    if (IGNORE_OMP || size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            size_t sz = mat.idi[i + 1] - mat.idi[i];
            auto matSpan = CSpan(mat.di, mat.idi[i], sz);
            SymmetricMatrix::multiplyMatrixBlockToVector(
                matSpan, x.cspan(i * 2, 2), res->span(i * 2, 2));

            for (size_t j = mat.ig[i]; j < mat.ig[i + 1]; ++j) {
                size_t k = mat.jg[j];
                size_t sz = mat.igg[j + 1] - mat.igg[j];

                auto matSpan = CSpan(mat.gg, mat.igg[j], sz);
                SymmetricMatrix::multiplyMatrixBlockToVector(
                    matSpan, x.cspan(k * 2, 2), res->span(i * 2, 2));
                SymmetricMatrix::multiplyMatrixBlockToVector(
                    matSpan, x.cspan(i * 2, 2), res->span(k * 2, 2));
            }
        }
    } else {
        std::vector<Vector> thread_results(NUM_THREADS);

#pragma omp parallel num_threads(NUM_THREADS)
        {
            auto t_num = omp_get_thread_num();
            auto& t_res = thread_results[t_num];
            t_res.resize(size * 2);

#pragma omp for schedule(dynamic, 300)
            for (size_t i = 0; i < size; ++i) {
                size_t sz = mat.idi[i + 1] - mat.idi[i];
                auto matSpan = CSpan(mat.di, mat.idi[i], sz);
                SymmetricMatrix::multiplyMatrixBlockToVector(
                    matSpan, x.cspan(i * 2, 2), t_res.span(i * 2, 2));

                for (size_t j = mat.ig[i]; j < mat.ig[i + 1]; ++j) {
                    size_t k = mat.jg[j];
                    size_t sz = mat.igg[j + 1] - mat.igg[j];
                    auto matSpan = CSpan(mat.gg, mat.igg[j], sz);
                    SymmetricMatrix::multiplyMatrixBlockToVector(
                        matSpan, x.cspan(k * 2, 2), t_res.span(i * 2, 2));
                    SymmetricMatrix::multiplyMatrixBlockToVector(
                        matSpan, x.cspan(i * 2, 2), t_res.span(k * 2, 2));
                }
            }
        }

#pragma omp parallel num_threads(NUM_THREADS)
        {
#pragma omp for VEC_SCHEDULE
            for (size_t j = 0; j < size * 2; ++j) {
                for (size_t i = 0; i < NUM_THREADS; ++i) {
                    auto& t_res = thread_results[i].data;
                    res->data[j] += t_res[j];
                }
            }
        }
    }
}

void SymmetricMatrix::getDiagonalPrecondition(const SymmetricMatrix& mat,
                                              Vector* diPrec) {
    auto size = mat.size();
    if (IGNORE_OMP || size < NUM_THREADS) {
        for (size_t i = 0; i < size; ++i) {
            int countItems = mat.idi[i + 1] - mat.idi[i];
            auto comp = std::complex<double>();

            comp.real(mat.di[mat.idi[i]]);

            if (countItems == 2) {
                comp.imag(mat.di[mat.idi[i]]);
            }

            comp = 1.0 / comp;

            diPrec->data[i * 2] = comp.real();
            diPrec->data[i * 2 + 1] = comp.imag();
        }
    } else {
#pragma omp parallel num_threads(NUM_THREADS)
        {
#pragma omp for VEC_SCHEDULE
            for (size_t i = 0; i < size; ++i) {
                int countItems = mat.idi[i + 1] - mat.idi[i];
                auto comp = std::complex<double>();

                comp.real(mat.di[mat.idi[i]]);

                if (countItems == 2) {
                    comp.imag(mat.di[mat.idi[i]]);
                }

                comp = 1.0 / comp;

                diPrec->data[i * 2] = comp.real();
                diPrec->data[i * 2 + 1] = comp.imag();
            }
        }
    }
}