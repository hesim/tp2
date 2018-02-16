cmake_minimum_required(VERSION 3.9)
project(tp2)

set(CMAKE_CXX_STANDARD 11)

add_executable(tp2
        PACC/Assert.hpp
        PACC/config.hpp
        PACC/Tokenizer.cpp
        PACC/Tokenizer.hpp
        Chrono.hpp
        lodepng.cpp
        lodepng.h
        tp2.cpp)
