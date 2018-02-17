cmake_minimum_required(VERSION 3.9)
project(tp2)

set(CMAKE_CXX_STANDARD 11)

add_custom_target(tp2 COMMAND make run -C ${tp2_SOURCE_DIR}
        CLION_EXE_DIR=${PROJECT_BINARY_DIR})