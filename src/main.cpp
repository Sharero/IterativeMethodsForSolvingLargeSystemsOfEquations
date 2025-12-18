#include <array>
#include <string>

#include "../include/slae.h"

int main() {
    for (const std::string& number :
         std::array<std::string, 3>{"1", "2", "3"}) {
        const std::string folder_index = number;

        SLAE slae;

        slae.inputComplexSLAEData(folder_index);
        slae.printComplexSLAEDataInformation();
        slae.solve(folder_index);
    }
}