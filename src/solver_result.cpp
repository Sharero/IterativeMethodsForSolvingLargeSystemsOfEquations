#include "../include/solver_result.h"

#include <iomanip>
#include <iostream>
#include <sstream>

std::string SolverResult::dump() const {
    std::ostringstream oss;

    oss << "{\n";

    oss << "  solved=" << (isSolved ? "true" : "false") << ",\n";

    oss << "  reachedEps=" << std::scientific << std::setprecision(14)
        << std::setw(20) << reachedEps << ",\n";

    oss << "  totalIterations=" << std::setw(8) << totalIterations << ",\n";

    oss << "  totalTimeSeconds=" << std::fixed << std::setprecision(14)
        << std::setw(20) << totalTimeSeconds << ",\n";

    oss << "}\n";

    return oss.str();
}

std::string SolverResult::dumpEps(std::string divider) const {
    std::ostringstream oss;

    oss << std::scientific << std::setprecision(14);

    for (size_t i = 0; i < eps.size(); ++i) {
        oss << std::setw(8) << i << divider << std::setw(20) << eps[i] << '\n';
    }

    return oss.str();
}