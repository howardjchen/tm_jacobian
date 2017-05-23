cmake_minimum_required(VERSION 2.8.2)
project(jacobian_gen)

find_package(Eigen3 REQUIRED)

add_executable(j jacobian_gen.cpp)
