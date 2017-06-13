cmake_minimum_required(VERSION 2.8.2)
project(jacobian_gen)

find_package(Eigen3 REQUIRED)

include_directories(tm_kinematics/include)

file(GLOB TM_KINEMATICS_SRC
	"tm_kinematics/src/tm_kin.cpp"
	"tm_kinematics/src/tm_jacob.cpp"
)

add_executable(tm_jacobian main.cpp ${TM_KINEMATICS_SRC})
