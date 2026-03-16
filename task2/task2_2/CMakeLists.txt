cmake_minimum_required(VERSION 3.10)
project(integral)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(integral integral.cpp)

target_link_libraries(integral PRIVATE OpenMP::OpenMP_CXX)