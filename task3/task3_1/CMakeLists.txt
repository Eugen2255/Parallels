cmake_minimum_required(VERSION 3.10)
project(dgemv CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3)

set(M "40000" CACHE STRING "Size of matrix dimension m")

add_executable(dgemv_thread dgemv_thread.cpp)
target_compile_definitions(dgemv_thread PRIVATE M=${M})

add_executable(dgemv_jthread dgemv_jthread.cpp)
target_compile_definitions(dgemv_jthread PRIVATE M=${M})

add_executable(dgemv_async dgemv_async.cpp)
target_compile_definitions(dgemv_async PRIVATE M=${M})


message(STATUS "Matrix size m = ${M}")
