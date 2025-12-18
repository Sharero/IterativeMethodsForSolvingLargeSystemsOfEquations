#ifndef COCG_DI_SOLVER_H
#define COCG_DI_SOLVER_H

#include "../include/solver_result.h"
#include "../include/symmetric_matrix.h"
#include "../include/vector.h"

class COCGDiSolver {
   public:
    SolverResult solve(const SymmetricMatrix& A, const Vector& b, Vector& x,
                       const bool isSmoothSolving);

    void setSolverParameters(int maximum_iterations_s, double epsilon_s,
                             int print_interval_s);

   private:
    int maximum_iterations = 0;
    double epsilon = 0.0;
    int print_interval = 0;

    SolverResult solveByBasicMethod(const SymmetricMatrix& A, const Vector& b,
                                    Vector& x);

    SolverResult solveBySmoothMethod(const SymmetricMatrix& A, const Vector& b,
                                     Vector& x);
};

#endif