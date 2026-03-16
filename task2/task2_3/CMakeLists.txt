cmake_minimum_required(VERSION 3.10)
project(equation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)

add_executable(equation_v1 var1.cpp)
add_executable(equation_v2 var2.cpp)

target_link_libraries(equation_v1 PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(equation_v2 PRIVATE OpenMP::OpenMP_CXX)