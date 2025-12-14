#include "../include/slae.h"

#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename T>
std::vector<T> ComplexSlae::readDataFromBinaryFile(
    const std::filesystem::path& file_name) {
    const auto file_size_u = std::filesystem::file_size(file_name);
    const auto count = static_cast<std::size_t>(file_size_u / sizeof(T));

    std::ifstream input_file(file_name, std::ios::binary);
    std::vector<char> buffer(static_cast<std::size_t>(file_size_u));

    input_file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));

    std::vector<T> result(count);

#pragma unroll 4
    for (std::size_t i = 0; i < count; ++i) {
        std::memcpy(&result.at(i), &buffer.at(i * sizeof(T)), sizeof(T));
    }

    return result;
}

void ComplexSlae::inputSLAEData(const std::string& folder_index) {
    const std::filesystem::path input_folder =
        std::filesystem::path("data") / "complex" / folder_index;

    const std::filesystem::path kuslau_path = input_folder / "kuslau";
    std::ifstream kuslau(kuslau_path);

    kuslau >> size >> epsilon >> maximum_iterations;
    kuslau.close();

    ig = readDataFromBinaryFile<int>(input_folder / "ig");
#pragma unroll 4
    for (auto& elem : ig) {
        elem -= 1;
    }
    std::vector<int> idi = readDataFromBinaryFile<int>(input_folder / "idi");

#pragma unroll 4
    for (auto& elem : idi) {
        elem -= 1;
    }

    const auto di_count = static_cast<std::size_t>(idi.at(size)) / 2;

    di = readDataFromBinaryFile<std::complex<double>>(input_folder / "di");
    jg = readDataFromBinaryFile<int>(input_folder / "jg");
#pragma unroll 4
    for (auto& elem : jg) {
        elem -= 1;
    }

    std::vector<int> ijg = readDataFromBinaryFile<int>(input_folder / "ijg");

#pragma unroll 4
    for (auto& elem : ijg) {
        elem -= 1;
    }

    const auto ig_size_index = static_cast<std::size_t>(ig.at(size));
    const auto gg_index =
        static_cast<std::size_t>(ijg.at(ig_size_index - 1)) / 2;

    gg = readDataFromBinaryFile<std::complex<double>>(input_folder / "gg");
    pr = readDataFromBinaryFile<std::complex<double>>(input_folder / "pr");

    x.resize(size, 0);

    std::ofstream out("data.txt");

#pragma unroll 4
    for (const auto elem : gg) {
        out << elem << "\n";
    }

    std::cout << "Read sizes:\n"
              << "  ig:  " << ig.size() << '\n'
              << "  idi: " << idi.size() << '\n'
              << "  di:  " << di.size() << "\n"
              << "  jg:  " << jg.size() << '\n'
              << "  ijg: " << ijg.size() << '\n'
              << "  gg:  " << gg.size() << "\n"
              << "  pr:  " << pr.size() << '\n'
              << "  x:  " << x.size() << "\n";

    ijg.clear();
    idi.clear();
}

void ComplexSlae::solveSLAE() {
    solver_COCG.init(ig, jg, di, gg, size);
    solver_COCG.solve(x, epsilon, pr, maximum_iterations);

    //     std::ofstream out("data.txt");

    // #pragma unroll 4
    //     for (const auto elem : x) {
    //         out << elem << "\n";
    //     }
}