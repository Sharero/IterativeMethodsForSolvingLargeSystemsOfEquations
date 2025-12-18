#ifndef SYMMETRIC_MATRIX_H
#define SYMMETRIC_MATRIX_H

#include <complex>
#include <optional>
#include <vector>

#include "../include/cspan.h"
#include "../include/span.h"
#include "../include/vector.h"

class SymmetricMatrix {
   public:
    std::vector<double> di;
    std::vector<double> gg;

    std::vector<int> idi;
    std::vector<int> igg;
    std::vector<int> ig;
    std::vector<int> jg;

    size_t size() const {
        return (idi.size() == 0) ? 0 : idi.size() - 1;
    }

    size_t totalItems() const {
        return (size() == 0) ? 0 : *ig.rbegin();
    }

    std::optional<std::complex<double>> tryGetEntry(size_t row,
                                                    size_t col) const;

    std::complex<double> at(size_t row, size_t col) const {
        auto val = tryGetEntry(row, col);
        return val.has_value() ? val.value() : std::complex<double>{};
    }

    static void multiplyMatrixBlockToVector(const CSpan& matSpan,
                                            const CSpan& vecSpan, Span resSpan);

    static void multiplyMatrixByVector(const SymmetricMatrix& mat,
                                       const Vector& x, Vector* res);

    static void getDiagonalPrecondition(const SymmetricMatrix& mat,
                                        Vector* diPrec);
};

#endif