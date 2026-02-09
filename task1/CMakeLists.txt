cmake_minimum_required(VERSION 3.10)
project(task1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_DOUBLE "double or float type" OFF)

if (USE_DOUBLE)
    add_definitions(-DDOUBLE)
    message(STATUS "Building with DOUBLE precision")
else()
    message(STATUS "Building with FLOAT precision") 
endif()

add_executable(task1 task1.cpp)

