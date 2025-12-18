#include "../include/COCG_Di_solver.h"

#include <iomanip>
#include <iostream>

#include "../include/solver_result.h"
#include "../include/timer.h"

void COCGDiSolver::setSolverParameters(int maximum_iterations_s,
                                       double epsilon_s, int print_interval_s) {
    maximum_iterations = maximum_iterations_s;
    print_interval = print_interval_s;
    epsilon = epsilon_s;
}

SolverResult COCGDiSolver::solve(const SymmetricMatrix& A, const Vector& b,
                                 Vector& x, const bool isSmoothSolving) {
    if (isSmoothSolving) {
        return solveBySmoothMethod(A, b, x);
    }

    return solveByBasicMethod(A, b, x);
}

SolverResult COCGDiSolver::solveByBasicMethod(const SymmetricMatrix& A,
                                              const Vector& b, Vector& x) {
    SolverResult _res;
    _res.eps.reserve(maximum_iterations);

    Vector r(b.size());
    Vector r0(b.size());
    Vector p(b.size());
    Vector di(b.size());
    Vector Ap(b.size());
    Vector z(b.size());

    std::complex<double> scal_r_z;
    double norm_r;
    double norm_b;

    Timer t;

    t.start();

    SymmetricMatrix::getDiagonalPrecondition(A, &di);
    SymmetricMatrix::multiplyMatrixByVector(A, x, &r);
    Vector::calculateLinearCombinationOfVectors(1.0, b, -1.0, r, &r);
    Vector::copyVector(r, &r0);
    Vector::multiplyVectorByVector(di, r, &z);
    Vector::copyVector(z, &p);

    scal_r_z = Vector::calculateVectorsProductWithoutConjugation(r, z);
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
        SymmetricMatrix::multiplyMatrixByVector(A, p, &Ap);

        auto alpha = Vector::calculateVectorsProductWithoutConjugation(r, z) /
                     Vector::calculateVectorsProductWithoutConjugation(Ap, p);

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, x, alpha, p,
                                                    &x);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, r, -alpha, Ap,
                                                    &r);
        auto norm_r = Vector::calculateNormOfRealPartOfVectorsProduct(r);
        eps = norm_r / norm_b;
        _res.eps.push_back(eps);

        Vector::multiplyVectorByVector(di, r, &z);

        auto next_scal_r_z =
            Vector::calculateVectorsProductWithoutConjugation(r, z);
        auto beta = next_scal_r_z / scal_r_z;
        scal_r_z = next_scal_r_z;

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, z, beta, p, &p);

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

SolverResult COCGDiSolver::solveBySmoothMethod(const SymmetricMatrix& A,
                                               const Vector& b, Vector& x) {
    SolverResult _res;
    _res.eps.reserve(maximum_iterations);

    Timer t;
    t.start();

    Vector r(b.size());
    Vector r0(b.size());
    Vector p(b.size());
    Vector di(b.size());
    Vector Ap(b.size());
    Vector z(b.size());
    Vector y(x.size());
    Vector s(r.size());
    Vector rs(r.size());

    std::complex<double> scal_r_z;
    double norm_s;
    double norm_b;

    SymmetricMatrix::getDiagonalPrecondition(A, &di);
    SymmetricMatrix::multiplyMatrixByVector(A, x, &r);
    Vector::calculateLinearCombinationOfVectors(1.0, b, -1.0, r, &r);
    Vector::copyVector(r, &r0);
    Vector::multiplyVectorByVector(di, r, &z);
    Vector::copyVector(z, &p);
    Vector::copyVector(x, &y);
    Vector::copyVector(r, &s);

    scal_r_z = Vector::calculateVectorsProductWithoutConjugation(r, z);
    norm_s = Vector::calculateNormOfRealPartOfVectorsProduct(s);
    norm_b = Vector::calculateNormOfRealPartOfVectorsProduct(b);

    double eps = norm_s / norm_b;
    size_t iteration = 0;

    std::cout << "Iter: " << std::setw(5) << iteration
              << ", eps: " << std::setprecision(14) << std::scientific << eps
              << '\n';

    _res.eps.push_back(eps);

    for (iteration = 1; iteration <= maximum_iterations && eps > epsilon;
         ++iteration) {
        SymmetricMatrix::multiplyMatrixByVector(A, p, &Ap);

        auto alpha = Vector::calculateVectorsProductWithoutConjugation(r, z) /
                     Vector::calculateVectorsProductWithoutConjugation(Ap, p);

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, x, alpha, p,
                                                    &x);
        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, r, -alpha, Ap,
                                                    &r);

        Vector::calculateLinearCombinationOfVectors(1.0, r, -1.0, s, &rs);
        double eta = -Vector::getRealPartOfVectorsProduct(s, rs) /
                     Vector::getRealPartOfVectorsProduct(rs, rs);
        eta = std::min(std::max(eta, 0.0), 1.0);
        Vector::calculateLinearCombinationOfVectors((1 - eta), y, eta, x, &y);
        Vector::calculateLinearCombinationOfVectors((1 - eta), s, eta, r, &s);

        auto norm_s = Vector::calculateNormOfRealPartOfVectorsProduct(s);
        eps = norm_s / norm_b;
        _res.eps.push_back(eps);

        Vector::multiplyVectorByVector(di, r, &z);

        auto next_scal_r_z =
            Vector::calculateVectorsProductWithoutConjugation(r, z);
        auto beta = next_scal_r_z / scal_r_z;
        scal_r_z = next_scal_r_z;

        Vector::calculateLinearCombinationOfVectors({1.0, 0.0}, z, beta, p, &p);

        if (iteration % print_interval == 0) {
            std::cout << "Iter: " << std::setw(5) << iteration
                      << ", eps: " << std::setprecision(14) << std::scientific
                      << eps << '\n';
        }
    }

    Vector::copyVector(y, &x);

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

    std::cout << "\n";

    return _res;
}