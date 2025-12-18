#ifndef VECTOR_H
#define VECTOR_H

#include <complex>
#include <vector>

#include "../include/cspan.h"
#include "../include/span.h"

class Vector {
   public:
    std::vector<double> data;

    Vector() = default;
    explicit Vector(size_t size) {
        resize(size);
    }

    void resize(size_t size) {
        data.resize(size);
    }
    size_t size() const {
        return data.size();
    }
    size_t blocks_count() const {
        return size() / 2;
    }

    std::complex<double> at(size_t ind) const {
        return {data.at(ind * 2), data.at(ind * 2 + 1)};
    }

    void setAt(size_t ind, const std::complex<double>& num) {
        data.at(ind * 2) = num.real();
        data.at(ind * 2 + 1) = num.imag();
    }

    void clear() {
        for (double& v : data) {
            v = 0.0;
        }
    }

    Span span(size_t start, size_t size) {
        return Span(data, start, size);
    }

    CSpan cspan(size_t start, size_t size) const {
        return CSpan(data, start, size);
    }

    static double getRealPartOfVectorsProduct(const Vector& x, const Vector& y);

    static std::complex<double> calculateVectorsProduct(const Vector& x,
                                                        const Vector& y);

    static std::complex<double> calculateVectorsProductWithoutConjugation(
        const Vector& x, const Vector& y);

    static double calculateNormOfRealPartOfVectorsProduct(const Vector& x);

    static void multiplyVectorByScalar(const Vector& x, double scal, Vector* y);

    static void multiplyVectorByVector(const Vector& x, const Vector& y,
                                       Vector* z);

    static void calculateLinearCombinationOfVectors(double mX, const Vector& x,
                                                    double mY, const Vector& y,
                                                    Vector* z);

    static void calculateLinearCombinationOfVectors(std::complex<double> mX,
                                                    const Vector& x,
                                                    std::complex<double> mY,
                                                    const Vector& y, Vector* z);

    static void copyVector(const Vector& x, Vector* y);
};

#endif