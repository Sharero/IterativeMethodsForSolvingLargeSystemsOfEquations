#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <complex>
#include <vector>

class MathUtils {
   public:
    static void multiplyMatrixByVector(
        const std::vector<std::complex<double>>& f,
        std::vector<std::complex<double>>& x,
        std::vector<std::complex<double>>& di, std::vector<int>& gi,
        std::vector<std::complex<double>>& gg, std::vector<int>& gj);

    static double calculateDotProduct(
        const std::vector<std::complex<double>>& vec);

    static std::complex<double> multiplyVectorByVector(
        const std::vector<std::complex<double>>& vec1,
        const std::vector<std::complex<double>>& vec2);
};

#endif