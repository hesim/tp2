cmake_minimum_required(VERSION 3.9)
project(tp2)

set(CMAKE_CXX_STANDARD 11)

include(FindOpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(tp2
        PACC/Assert.hpp
        PACC/config.hpp
        PACC/Tokenizer.cpp
        PACC/Tokenizer.hpp
        Chrono.hpp
        lodepng.cpp
        lodepng.h
        tp2.cpp)
