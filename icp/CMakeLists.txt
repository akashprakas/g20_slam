cmake_minimum_required(VERSION 3.22)
project(icp)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_library(icp icp2d.cpp icp2d.h)

add_executable(test test.cpp )

target_link_libraries(test icp)