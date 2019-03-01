cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(MAGPROF)
enable_language(Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
#include_directories(${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(GLOB SRC_FILES "src/*.f90")

add_executable(MAGPROF ${SRC_FILES})
