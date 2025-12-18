#include "../include/vector.h"

#include <omp.h>

static constexpr bool IGNORE_OMP = false;
static constexpr size_t NUM_THREADS = 8;
static constexpr size_t MIN_SIZE = 1000;

#define VEC_SCHEDULE schedule(static)

double Vector::getRealPartOfVectorsProduct(const Vector& x, const Vector& y) {
    return Vector::calculateVectorsProduct(x, y).real();
}

std::complex<double> Vector::calculateVectorsProduct(const Vector& x,
                                                     const Vector& y) {
    size_t size = x.blocks_count();
    std::complex<double> res;

    if (IGNORE_OMP || size * size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            std::complex<double> comp_l = {
                x.data[i * 2],
                -x.data[i * 2 + 1],
            };

            std::complex<double> comp_r = {
                y.data[i * 2],
                y.data[i * 2 + 1],
            };

            res += comp_l * comp_r;
        }

        return res;
    }
    double real_acc = 0.0;
    double imag_acc = 0.0;

#pragma omp parallel num_threads(NUM_THREADS)
    {
#pragma omp for reduction(+ : real_acc, imag_acc) VEC_SCHEDULE
        for (size_t i = 0; i < size; ++i) {
            const double xr = x.data[2 * i];
            const double xi = x.data[2 * i + 1];
            const double yr = y.data[2 * i];
            const double yi = y.data[2 * i + 1];

            real_acc += xr * yr + xi * yi;
            imag_acc += xr * yi - xi * yr;
        }
    }

    return std::complex<double>(real_acc, imag_acc);
}

std::complex<double> Vector::calculateVectorsProductWithoutConjugation(
    const Vector& x, const Vector& y) {
    auto size = x.blocks_count();
    std::complex<double> res;

    if (IGNORE_OMP || size * size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            std::complex<double> comp_l = {
                x.data[i * 2],
                x.data[i * 2 + 1],
            };

            std::complex<double> comp_r = {
                y.data[i * 2],
                y.data[i * 2 + 1],
            };

            res += comp_l * comp_r;
        }

        return res;
    }

    double real_acc = 0.0;
    double imag_acc = 0.0;

#pragma omp parallel num_threads(NUM_THREADS)
    {
#pragma omp for reduction(+ : real_acc, imag_acc) VEC_SCHEDULE
        for (size_t i = 0; i < size; ++i) {
            const double xr = x.data[2 * i];
            const double xi = x.data[2 * i + 1];
            const double yr = y.data[2 * i];
            const double yi = y.data[2 * i + 1];

            real_acc += xr * yr - xi * yi;
            imag_acc += xr * yi + xi * yr;
        }
    }

    return std::complex<double>(real_acc, imag_acc);
}

double Vector::calculateNormOfRealPartOfVectorsProduct(const Vector& x) {
    return std::sqrt(getRealPartOfVectorsProduct(x, x));
}

void Vector::multiplyVectorByScalar(const Vector& x, double scal, Vector* y) {
    size_t size = x.size();
    std::vector<double>& yy = y->data;
    const std::vector<double>& xx = x.data;

    if (IGNORE_OMP || size * size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            yy[i] = xx[i] * scal;
        }
    } else {
#pragma omp parallel num_threads(NUM_THREADS)
        {
#pragma omp for VEC_SCHEDULE
            for (size_t i = 0; i < size; ++i) {
                yy[i] = xx[i] * scal;
            }
        }
    }
}

void Vector::multiplyVectorByVector(const Vector& x, const Vector& y,
                                    Vector* z) {
    size_t size = x.blocks_count();

    if (IGNORE_OMP || size * size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            std::complex<double> xx = {x.data[i * 2], x.data[i * 2 + 1]};
            std::complex<double> yy = {y.data[i * 2], y.data[i * 2 + 1]};
            std::complex<double> res = xx * yy;

            z->data[i * 2] = res.real();
            z->data[i * 2 + 1] = res.imag();
        }
    } else {
#pragma omp parallel num_threads(NUM_THREADS)
        {
#pragma omp for VEC_SCHEDULE
            for (size_t i = 0; i < size; ++i) {
                std::complex<double> xx = {x.data[i * 2], x.data[i * 2 + 1]};
                std::complex<double> yy = {y.data[i * 2], y.data[i * 2 + 1]};
                std::complex<double> res = xx * yy;

                z->data[i * 2] = res.real();
                z->data[i * 2 + 1] = res.imag();
            }
        }
    }
}

void Vector::calculateLinearCombinationOfVectors(double mX, const Vector& x,
                                                 double mY, const Vector& y,
                                                 Vector* z) {
    size_t size = x.size();
    const std::vector<double>& xx = x.data;
    const std::vector<double>& yy = y.data;
    std::vector<double>& zz = z->data;

    if (IGNORE_OMP || size * size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            zz[i] = mX * xx[i] + mY * yy[i];
        }
    } else {
#pragma omp parallel num_threads(NUM_THREADS)
        {
#pragma omp for VEC_SCHEDULE
            for (size_t i = 0; i < size; ++i) {
                zz[i] = mX * xx[i] + mY * yy[i];
            }
        }
    }
}

void Vector::calculateLinearCombinationOfVectors(std::complex<double> mX,
                                                 const Vector& x,
                                                 std::complex<double> mY,
                                                 const Vector& y, Vector* z) {
    size_t size = x.blocks_count();
    const std::vector<double>& xx = x.data;
    const std::vector<double>& yy = y.data;
    std::vector<double>& zz = z->data;

    if (IGNORE_OMP || size * size < MIN_SIZE) {
        for (size_t i = 0; i < size; ++i) {
            std::complex<double> x_comp = {xx[i * 2], xx[i * 2 + 1]};
            std::complex<double> y_comp = {yy[i * 2], yy[i * 2 + 1]};
            std::complex<double> res = mX * x_comp + mY * y_comp;

            zz[i * 2] = res.real();
            zz[i * 2 + 1] = res.imag();
        }
    } else {
#pragma omp parallel num_threads(NUM_THREADS)
        {
#pragma omp for VEC_SCHEDULE
            for (size_t i = 0; i < size; ++i) {
                std::complex<double> x_comp = {xx[i * 2], xx[i * 2 + 1]};
                std::complex<double> y_comp = {yy[i * 2], yy[i * 2 + 1]};
                std::complex<double> res = mX * x_comp + mY * y_comp;

                zz[i * 2] = res.real();
                zz[i * 2 + 1] = res.imag();
            }
        }
    }
}

void Vector::copyVector(const Vector& x, Vector* y) {
    Vector::multiplyVectorByScalar(x, 1.0, y);
}