#ifndef SLAE_H
#define SLAE_H

#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "COCG.h"

class ComplexSlae {
   private:
    int size = 0;
    int maximum_iterations = 0;
    double epsilon = 0.0;

    std::vector<int> ig;
    std::vector<int> idi;
    std::vector<int> jg;
    std::vector<int> ijg;

    std::vector<std::complex<double>> di;
    std::vector<std::complex<double>> gg;
    std::vector<std::complex<double>> pr;
    std::vector<std::complex<double>> x;

    SolverCOCG solver_COCG;

    template <typename T>
    std::vector<T> readDataFromBinaryFile(
        const std::filesystem::path& file_name);

   public:
    void inputSLAEData(const std::string& folder_index);
    void solveSLAE();
};

#endif