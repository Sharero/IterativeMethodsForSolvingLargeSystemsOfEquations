cmake_minimum_required(VERSION 3.20)

project(IterativeMethodsForSolvingLargeSystemsOfEquations
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

add_executable(
    my_app
    src/main.cpp
    src/slae.cpp
    src/vector.cpp
    src/symmetric_matrix.cpp
    src/COCG_Di_solver.cpp
    src/COCR_Di_solver.cpp
    src/solver_result.cpp
)

target_link_libraries(my_app PRIVATE OpenMP::OpenMP_CXX)
