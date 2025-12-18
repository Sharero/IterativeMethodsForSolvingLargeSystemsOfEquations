#include "../include/COCR_Di_solver.h"

#include <iomanip>
#include <iostream>

#include "../include/solver_result.h"
#include "../include/timer.h"

void COCRDiSolver::setSolverParameters(int maximum_iterations_s,
                                       double epsilon_s, int print_interval_s) {
    maximum_iterations = maximum_iterations_s;
    print_interval = print_interval_s;
    epsilon = epsilon_s;
}

SolverResult COCRDiSolver::solve(const SymmetricMatrix& A, const Vector& b,
                                 Vector& x, const bool isSmoothSolving) {
    if (isSmoothSolving) {
        return solveBySmoothMethod(A, b, x);
    }

    return solveByBasicMethod(A, b, x);
}

SolverResult COCRDiSolver::solveByBasicMethod(const SymmetricMatrix& A,
                                              const Vector& b, Vector& x) {
    SolverResult _res;
    _res.eps.reserve(maximum_iterations);

    Vector r(b.size());
    Vector r0(b.size());
    Vector p(b.size());
    Vector di(b.size());
    Vector Ap(b.size());
    Vector z(b.size());
    Vector s(b.size());
    Vector a(b.size());
    Vector w(b.size());

    std::complex<double> scal_a_s;
    double norm_r;
    double norm_b;

    Timer t;

    t.start();

    SymmetricMatrix::getDiagonalPrecondition(A, &di);  // di = 1.0 / A.di

    SymmetricMatrix::multiplyMatrixByVector(A, x, &r);  // r = b - A*x
    Vector::calculateLinearCombinationOfVectors(1.0, b, -1.0, r,
                                                &r);  // r = b - A*x

    Vector::copyVector(r, &r0);                 // r0 = r
    Vector::multiplyVectorByVector(di, r, &s);  // s = di^-1 * r
    Vector::copyVector(s, &p);                  // p = s

    SymmetricMatrix::multiplyMatrixByVector(A, s, &z);  // z = A*s
    Vector::copyVector(z, &a);                          // a = z

    Vector::multiplyVectorByVector(di, z, &w);  // w = di^-1 * z

    scal_a_s = Vector::calculateVectorsProductWithoutConjugation(a, s);
    norm_r = Vector::calculateNormOfRealPartOfVectorsProduct(r);
    norm_b = Vector::calculateNormOfRealPartOfVectorsProduct(b);

    double eps = norm_r / norm_b;
    size_t iteration = 0;

    std::cout << "Iter: " << std::setw(5) << iteration
              << ", eps: " << std::setprecision(14) << std::scientific << eps
              << '\n';

    _res.eps.push_back(eps);

    for (iteration = 1; iteration <= maximum_iterations && eps > epsilon;
         ++iteration) {
        auto alpha = Vector::calculateVectorsProductWithoutConjugation(a, s) /
                     Vector::calculateVectorsProductWithoutConjugation(z, w);

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, x, alpha, p,
                                                    &x);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, r, -alpha, z,
                                                    &r);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, s, -alpha, w,
                                                    &s);

        SymmetricMatrix::multiplyMatrixByVector(A, s, &a);

        auto norm_r = Vector::calculateNormOfRealPartOfVectorsProduct(r);
        eps = norm_r / norm_b;
        _res.eps.push_back(eps);

        auto next_scal_a_s =
            Vector::calculateVectorsProductWithoutConjugation(a, s);

        auto beta = next_scal_a_s / scal_a_s;

        scal_a_s = next_scal_a_s;

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, s, beta, p, &p);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, a, beta, z, &z);

        Vector::multiplyVectorByVector(di, z, &w);

        if (iteration % print_interval == 0) {
            std::cout << "Iter: " << std::setw(5) << iteration
                      << ", eps: " << std::setprecision(14) << std::scientific
                      << eps << '\n';
        }
    }

    t.stop();

    _res.reachedEps = eps;
    _res.totalIterations = iteration - 1;
    _res.totalTimeSeconds = t.elapsedSeconds();
    _res.isSolved = true;

    std::cout << "COCG diag precond solver finished\n";
    std::cout << " - Elapsed time: " << t.elapsedSeconds() << " seconds\n";
    std::cout << " - Total iterations: " << std::setw(6) << (iteration - 1)
              << " (limit: " << std::setw(6) << maximum_iterations << ")\n";
    std::cout << " - Reached epsilon: " << std::setw(6) << std::scientific
              << std::setprecision(2) << eps << " (target: " << std::setw(6)
              << epsilon << ")\n";

    SymmetricMatrix::multiplyMatrixByVector(A, x, &r);
    Vector::calculateLinearCombinationOfVectors(1.0, b, -1.0, r, &r);
    auto norm_residual = Vector::calculateNormOfRealPartOfVectorsProduct(r) /
                         Vector::calculateNormOfRealPartOfVectorsProduct(r0);
    std::cout << " - Norm of relative residual: " << std::scientific
              << std::setprecision(2) << norm_residual << "\n";

    if (iteration - 1 == maximum_iterations && std::abs(eps / epsilon) > 5) {
        std::cout << "Warning: Solver didn't reach target epsilon.\n";
        _res.isSolved = false;
    }

    return _res;
}

SolverResult COCRDiSolver::solveBySmoothMethod(const SymmetricMatrix& A,
                                               const Vector& b, Vector& x) {
    SolverResult _res;
    _res.eps.reserve(maximum_iterations);

    Vector r(b.size());
    Vector r0(b.size());
    Vector p(b.size());
    Vector di(b.size());
    Vector Ap(b.size());
    Vector z(b.size());
    Vector s(b.size());
    Vector a(b.size());
    Vector w(b.size());

    Vector resid_y(b.size());
    Vector resid_s(b.size());

    Vector r_resid_s(b.size());

    std::complex<double> scal_a_s;
    double norm_r;
    double norm_b;

    Timer t;

    t.start();

    SymmetricMatrix::getDiagonalPrecondition(A, &di);  // di = 1.0 / A.di

    SymmetricMatrix::multiplyMatrixByVector(A, x, &r);  // r = b - A*x
    Vector::calculateLinearCombinationOfVectors(1.0, b, -1.0, r,
                                                &r);  // r = b - A*x

    Vector::copyVector(r, &r0);                 // r0 = r
    Vector::multiplyVectorByVector(di, r, &s);  // s = di^-1 * r
    Vector::copyVector(s, &p);                  // p = s

    SymmetricMatrix::multiplyMatrixByVector(A, s, &z);  // z = A*s
    Vector::copyVector(z, &a);                          // a = z

    Vector::multiplyVectorByVector(di, z, &w);  // w = di^-1 * z

    Vector::copyVector(x, &resid_y);
    Vector::copyVector(r, &resid_s);

    scal_a_s = Vector::calculateVectorsProductWithoutConjugation(a, s);
    norm_r = Vector::calculateNormOfRealPartOfVectorsProduct(resid_s);
    norm_b = Vector::calculateNormOfRealPartOfVectorsProduct(b);

    double eps = norm_r / norm_b;
    size_t iteration = 0;

    std::cout << "Iter: " << std::setw(5) << iteration
              << ", eps: " << std::setprecision(14) << std::scientific << eps
              << '\n';

    _res.eps.push_back(eps);

    for (iteration = 1; iteration <= maximum_iterations && eps > epsilon;
         ++iteration) {
        auto alpha =
            scal_a_s / Vector::calculateVectorsProductWithoutConjugation(z, w);

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, x, alpha, p,
                                                    &x);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, r, -alpha, z,
                                                    &r);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, s, -alpha, w,
                                                    &s);

        Vector::calculateLinearCombinationOfVectors(1.0, r, -1.0, resid_s,
                                                    &r_resid_s);
        double eta = -Vector::getRealPartOfVectorsProduct(resid_s, r_resid_s) /
                     Vector::getRealPartOfVectorsProduct(r_resid_s, r_resid_s);
        eta = std::min(std::max(eta, 0.0), 1.0);
        Vector::calculateLinearCombinationOfVectors((1 - eta), resid_y, eta, x,
                                                    &resid_y);
        Vector::calculateLinearCombinationOfVectors((1 - eta), resid_s, eta, r,
                                                    &resid_s);

        SymmetricMatrix::multiplyMatrixByVector(A, s, &a);

        auto norm_s = Vector::calculateNormOfRealPartOfVectorsProduct(resid_s);
        eps = norm_s / norm_b;
        _res.eps.push_back(eps);

        auto next_scal_a_s =
            Vector::calculateVectorsProductWithoutConjugation(a, s);

        auto beta = next_scal_a_s / scal_a_s;

        scal_a_s = next_scal_a_s;

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, s, beta, p, &p);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, a, beta, z, &z);

        Vector::multiplyVectorByVector(di, z, &w);

        if (iteration % print_interval == 0) {
            std::cout << "Iter: " << std::setw(5) << iteration
                      << ", eps: " << std::setprecision(14) << std::scientific
                      << eps << '\n';
        }
    }

    t.stop();

    Vector::copyVector(resid_y, &x);

    _res.reachedEps = eps;
    _res.totalIterations = iteration - 1;
    _res.totalTimeSeconds = t.elapsedSeconds();
    _res.isSolved = true;

    std::cout << "COCG diag precond solver finished\n";
    std::cout << " - Elapsed time: " << t.elapsedSeconds() << " seconds\n";
    std::cout << " - Total iterations: " << std::setw(6) << (iteration - 1)
              << " (limit: " << std::setw(6) << maximum_iterations << ")\n";
    std::cout << " - Reached epsilon: " << std::setw(6) << std::scientific
              << std::setprecision(2) << eps << " (target: " << std::setw(6)
              << epsilon << ")\n";

    SymmetricMatrix::multiplyMatrixByVector(A, x, &r);
    Vector::calculateLinearCombinationOfVectors(1.0, b, -1.0, r, &r);
    auto norm_residual = Vector::calculateNormOfRealPartOfVectorsProduct(r) /
                         Vector::calculateNormOfRealPartOfVectorsProduct(r0);
    std::cout << " - Norm of relative residual: " << std::scientific
              << std::setprecision(2) << norm_residual << "\n";

    if (iteration - 1 == maximum_iterations && std::abs(eps / epsilon) > 5) {
        std::cout << "Warning: Solver didn't reach target epsilon.\n";
        _res.isSolved = false;
    }

    return _res;
}