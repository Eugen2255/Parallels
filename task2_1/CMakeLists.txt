cmake_minimum_required(VERSION 3.10)
project(dgemv CXX)

set(M "40000" CACHE STRING "Size of matrix dimension m")

find_package(OpenMP REQUIRED)

add_executable(dgemv dgemv.cpp)
target_compile_definitions(dgemv PRIVATE M=${M})
target_link_libraries(dgemv PRIVATE OpenMP::OpenMP_CXX)

set_target_properties(dgemv PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

#target_compile_options(dgemv PRIVATE 
#    -O1  
#)

message(STATUS "Matrix size m = ${M}")
