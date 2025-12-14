#include "../include/COCG.h"

#include <complex>
#include <iostream>
#include <vector>

#include "../include/mathUtils.h"

void SolverCOCG::init(std::vector<int>& gi_s, std::vector<int>& gj_s,
                      std::vector<std::complex<double>>& di_s,
                      std::vector<std::complex<double>>& gg_s, int n_s) {
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    r.resize(n);
    z.resize(n);
    p.resize(n);
    s.resize(n);
}

void SolverCOCG::solve(std::vector<std::complex<double>>& solution, double eps,
                       std::vector<std::complex<double>>& rp_s,
                       int maximum_iterations) {
    const int print_interval = 10;

    rp = rp_s;

    MathUtils::multiplyMatrixByVector(solution, r, di, gi, gg, gj);

//     x0.resize(n);

// #pragma unroll 4
//     for (int i = 0; i < n; i++) {
//         x0[i] = solution[i];
//         r[i] = rp[i] - r[i];
//         z[i] = r[i];
//         p[i] = z[i];
//     }

//     std::complex<double> alpha;
//     std::complex<double> beta;
//     std::complex<double> prod_1;
//     std::complex<double> prod_2;

//     double dicrepancy = 0.0;
//     double rp_norm = sqrt(MathUtils::calculateDotProduct(rp));

//     if (std::isnan(rp_norm)) {
// #pragma unroll 4
//         for (int i = 0; i < n; i++) {
//             solution[i] = x0[i];
//         }

//         x0.clear();
//         return;
//     }

//     prod_1 = MathUtils::multiplyVectorByVector(p, r);

//     bool finished = false;

//     int iteration = 0;
//     for (iteration = 0; iteration <= maximum_iterations && !finished;
//          iteration++) {
//         dicrepancy = sqrt(MathUtils::calculateDotProduct(r));

//         if (std::isnan(dicrepancy)) {
// #pragma unroll 4
//             for (int i = 0; i < n; i++) {
//                 solution[i] = x0[i];
//             }

//             x0.clear();
//             return;
//         }

//         double residual = dicrepancy / rp_norm;

//         if (iteration % print_interval == 0) {
//             std::cout << "COCG Residual: " << iteration << " " << residual
//                       << "\n";
//         }

//         if (residual > eps) {
//             MathUtils::multiplyMatrixByVector(z, s, di, gi, gg, gj);

//             alpha = prod_1 / MathUtils::multiplyVectorByVector(s, z);

// #pragma unroll 4
//             for (int i = 0; i < n; i++) {
//                 x0[i] += alpha * z[i];
//                 r[i] -= alpha * s[i];
//                 p[i] = r[i];
//             }

//             prod_2 = MathUtils::multiplyVectorByVector(p, r);

//             beta = prod_2 / prod_1;

//             prod_1 = prod_2;

// #pragma unroll 4
//             for (int i = 0; i < n; i++) {
//                 z[i] = p[i] + beta * z[i];
//             }
//         } else {
//             finished = true;
//         }
//     }

//     MathUtils::multiplyMatrixByVector(x0, r, di, gi, gg, gj);

// #pragma unroll 4
//     for (size_t i = 0; i < n; i++) {
//         r[i] = rp[i] - r[i];
//     }
//     dicrepancy = sqrt(MathUtils::calculateDotProduct(r));

//     std::cout << "COCG Residual: " << iteration - 1 << " "
//               << dicrepancy / rp_norm << "\n";

//     if (iteration >= maximum_iterations) {
//         std::cout << "Soulution can`t found, iteration limit exceeded!" << "\n";
//     }

// #pragma unroll 4
//     for (int i = 0; i < n; i++) {
//         solution[i] = x0[i];
//     }

//     x0.clear();
}
