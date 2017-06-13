cmake_minimum_required(VERSION 2.8.2)
project(jacobian_gen)

find_package(Eigen3 REQUIRED)

include_directories(tm_kinematics)

file(GLOB TM_KINEMATICS_LIB
	"tm_kinematics/tm_kinematics/tm_kin.h"
	"tm_kinematics/tm_kin.cpp"
	"tm_kinematics/tm_kinematics/tm_jacob.h"
	"tm_kinematics/tm_jacob.cpp"
)

add_executable(tm_jacobian main.cpp ${TM_KINEMATICS_LIB})
