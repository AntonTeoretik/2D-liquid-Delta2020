project(project)
cmake_minimum_required(VERSION 2.8)
set(SOURCES main.cpp bitmap_image.hpp diffeq_solvers.cpp diffeq_solvers.h diff_operators.cpp diff_operators.h drawer.cpp drawer.h grids.cpp grids.h logger.cpp logger.h matrix.cpp matrix.h slae_solvers.cpp slae_solvers.h)
set(CMAKE_CXX_FLAGS -std=c++11)
aux_source_directory(. SRC_LIST)
include_directories(.)
add_executable(${PROJECT_NAME} ${SOURCES})

