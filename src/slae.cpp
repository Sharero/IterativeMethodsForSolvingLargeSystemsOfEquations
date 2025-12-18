#include "../include/slae.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

void SLAE::solve(const std::string& folder_index) {
    COCG_Di_solver.setSolverParameters(maximum_iterations, epsilon, 10);
    COCR_Di_solver.setSolverParameters(maximum_iterations, epsilon, 10);

    // auto start_x_for_COCG = x;
    // auto COCG_results = COCG_Di_solver.solve(matrix, b, start_x_for_COCG,
    // false); auto end_x_for_COCG = x;
    // SymmetricMatrix::multiplyMatrixByVector(matrix, start_x_for_COCG,
    //                                         &end_x_for_COCG);

    // auto start_x_for_COCG_smooth = x;
    // auto COCG_smooth_results =
    //     COCG_Di_solver.solve(matrix, b, start_x_for_COCG_smooth, true);
    // auto end_x_for_COCG_smooth = x;
    // SymmetricMatrix::multiplyMatrixByVector(matrix, start_x_for_COCG_smooth,
    //                                         &end_x_for_COCG_smooth);

    // auto start_x_for_COCR = x;
    // auto COCR_results = COCR_Di_solver.solve(matrix, b, start_x_for_COCR,
    // false); auto end_x_for_COCR = x;
    // SymmetricMatrix::multiplyMatrixByVector(matrix, start_x_for_COCR,
    //                                         &end_x_for_COCR);

    auto start_x_for_COCR_smooth = x;
    auto COCR_smooth_results =
        COCR_Di_solver.solve(matrix, b, start_x_for_COCR_smooth, true);
    auto end_x_for_COCR_smooth = x;
    SymmetricMatrix::multiplyMatrixByVector(matrix, start_x_for_COCR_smooth,
                                            &end_x_for_COCR_smooth);

    // {
    //   auto out = std::ofstream("output/omp8/COCG_Di_" + folder_index +
    //   ".txt"); out << COCG_results.dump() << "\n"; out <<
    //   COCG_results.dumpEps() << "\n";
    // }

    // {
    //   auto out =
    //       std::ofstream("output/omp8/COCG_Di_Smooth_" + folder_index +
    //       ".txt");
    //   out << COCG_smooth_results.dump() << "\n";
    //   out << COCG_smooth_results.dumpEps() << "\n";
    // }

    // {
    //   auto out = std::ofstream("output/omp8/COCR_Di_" + folder_index +
    //   ".txt"); out << COCR_results.dump() << "\n"; out <<
    //   COCR_results.dumpEps() << "\n";
    // }

    {
        auto out = std::ofstream("output/omp8/COCR_Di_Smooth_" + folder_index +
                                 ".txt");
        out << COCR_smooth_results.dump() << "\n";
        out << COCR_smooth_results.dumpEps() << "\n";
    }
}

template <typename T>
void SLAE::readBinaryFileOfData(const std::filesystem::path& file_name,
                                T* result, int number_of_records,
                                int len_of_record) {
    std::ifstream in(file_name, std::ios::binary);

    in.read(reinterpret_cast<char*>(result),
            static_cast<std::streamsize>(sizeof(T) * len_of_record *
                                         number_of_records));
}

void SLAE::printComplexSLAEDataInformation() {
    std::cout << "SLAE information:" << "\n";
    std::cout << "\tdi: " << matrix.di.size() << "\n";
    std::cout << "\tgg: " << matrix.gg.size() << "\n";
    std::cout << "\tidi: " << matrix.idi.size() << "\n";
    std::cout << "\tig: " << matrix.ig.size() << "\n";
    std::cout << "\tigg: " << matrix.igg.size() << "\n";
    std::cout << "\tjg: " << matrix.jg.size() << "\n";
    std::cout << "\tx: " << x.size() << "\n";
    std::cout << "\tb: " << b.size() << "\n";
}

void SLAE::inputComplexSLAEData(const std::string& folder_index) {
    const std::filesystem::path input_folder =
        std::filesystem::path("data") / "complex" / folder_index;

    std::vector<int> ig;
    std::vector<int> igg;
    std::vector<int> jg;
    std::vector<int> idi;
    std::vector<double> di;
    std::vector<double> gg;
    std::vector<double> pr;

    {
        std::ifstream kuslau(input_folder / "kuslau");
        kuslau >> size >> epsilon >> maximum_iterations;
    }

    ig.resize(size + 1);
    readBinaryFileOfData(input_folder / "ig", ig.data(), ig.size(), 1);

    for (auto& i : ig) {
        --i;
    }

    idi.resize(size + 1);
    readBinaryFileOfData(input_folder / "idi", idi.data(), idi.size(), 1);

    for (auto& i : idi) {
        --i;
    }

    int ig_n_1 = ig[size];
    int di_n_1 = idi[size];

    jg.resize(ig_n_1);
    readBinaryFileOfData(input_folder / "jg", jg.data(), jg.size(), 1);

    for (auto& j : jg) {
        --j;
    }

    igg.resize(ig_n_1 + 1);
    readBinaryFileOfData(input_folder / "ijg", igg.data(), igg.size(), 1);

    for (auto& i : igg) {
        --i;
    }

    di.resize(di_n_1);
    readBinaryFileOfData(input_folder / "di", di.data(), di.size(), 1);

    int gg_count = igg[ig_n_1];

    gg.resize(gg_count);
    readBinaryFileOfData(input_folder / "gg", gg.data(), gg.size(), 1);

    pr.resize(size * 2);
    readBinaryFileOfData(input_folder / "pr", pr.data(), pr.size(), 1);

    matrix.di = std::move(di);
    matrix.gg = std::move(gg);
    matrix.idi = std::move(idi);
    matrix.ig = std::move(ig);
    matrix.jg = std::move(jg);
    matrix.igg = std::move(igg);

    b.data = std::move(pr);
    x.resize(b.size());
}