cmake_minimum_required(VERSION 2.8.2)
project(jacobian_gen)

include_directories(tm_kinematics)

file(GLOB TM_KINEMATICS_LIB
	"tm_kinematics/tm_kinematics/tm_kin.h"
	"tm_kinematics/tm_kin.cpp"
)

find_package(Eigen3 REQUIRED)

add_executable(tm_jacobian jacobian_gen.cpp ${TM_KINEMATICS_LIB})
