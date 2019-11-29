cmake_minimum_required(VERSION 3.14)
project(JacobiSVD)

set(CMAKE_CXX_STANDARD 14)

add_executable(JacobiSVD main.cpp Matrix/Matrix.cpp Matrix/Matrix.h common.h Serial/Serial.cpp Serial/Serial.h Parallel/Parallel.cpp Parallel/Parallel.h BlockParallel/BlockParallel.cpp BlockParallel/BlockParallel.h)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()