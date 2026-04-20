cmake_minimum_required(VERSION 3.10)
project(main CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#add_compile_options(-O2)

add_executable(server main.cpp)