#ifndef SOLVER_RESULT_H
#define SOLVER_RESULT_H

#include <string>
#include <vector>

class SolverResult {
   public:
    bool isSolved = false;
    double reachedEps = 0.0;
    size_t totalIterations = 0;
    double totalTimeSeconds = 0;

    std::vector<double> eps;

    std::string dump() const;
    std::string dumpEps(std::string divider = " ; ") const;
};

#endif