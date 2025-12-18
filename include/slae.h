#ifndef SLAE_H
#define SLAE_H

#include <filesystem>
#include <string>

#include "../include/COCG_Di_solver.h"
#include "../include/COCR_Di_solver.h"
#include "../include/symmetric_matrix.h"
#include "../include/vector.h"

class SLAE {
   private:
    int size = 0;
    double epsilon = 0.0;
    int maximum_iterations = 0;

    SymmetricMatrix matrix;
    Vector x;
    Vector b;

    COCGDiSolver COCG_Di_solver;
    COCRDiSolver COCR_Di_solver;

    template <typename T>
    void readBinaryFileOfData(const std::filesystem::path& file_name, T* result,
                              int number_of_records, int len_of_record);

   public:
    void solve(const std::string& folder_index);

    void inputComplexSLAEData(const std::string& folder_index);

    void printComplexSLAEDataInformation();
};

#endif