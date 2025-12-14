#include "../include/mathUtils.h"

#include <iostream>
void MathUtils::multiplyMatrixByVector(
    const std::vector<std::complex<double>>& f,
    std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& di,
    std::vector<int>& gi, std::vector<std::complex<double>>& gg,
    std::vector<int>& gj) {
    for (int i = 0; i < x.size(); i++) {
        std::complex<double> temp_elem = f[i];
        x[i] = di[i] * temp_elem;

#pragma unroll 4
        for (int k = gi[i], k1 = gi[i + 1]; k < k1; k++) {
            int j = gj[k];
            std::cout << i << " " << k << " " << gg[k] << "\n";
            // std::cout << i << " " << x[i] << " " << gg[k] << " " << f[j]
            //           << "\n";
            // x[i] += gg[k] * f[j];
            // x[j] += gg[k] * temp_elem;
        }
    }
}

double MathUtils::calculateDotProduct(
    const std::vector<std::complex<double>>& vec) {
    double result = 0.0;

#pragma unroll 4
    for (auto i : vec) {
        double re = i.real();
        double im = i.imag();
        result += (re * re) + (im * im);
    }

    return result;
}

std::complex<double> MathUtils::multiplyVectorByVector(
    const std::vector<std::complex<double>>& vec1,
    const std::vector<std::complex<double>>& vec2) {
    std::complex<double> result = 0.0;

#pragma unroll 4
    for (int i = 0; i < vec1.size(); i++) {
        result += vec1[i] * vec2[i];
    }

    return result;
}