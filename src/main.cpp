#include <string>

#include "../include/slae.h"

int main() {
    const std::string folder_index = "1";

    ComplexSlae complex_slae;

    complex_slae.inputSLAEData(folder_index);
}