/*********************************************************************
 * jacobian_cpp.cpp
 *
 * Copyright (c) 2017, ISCI / National Chiao Tung University (NCTU)
 *
 * Author: Howard Chen (s880367@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **********************************************************************/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include "tm_kinematics/include/tm_kinematics/tm_kin.h"
#include <unistd.h>
#include <sys/time.h>


#define D1 0.1451
#define D2 0
#define D3 0
#define D4 -0.1222
#define D5 0.106
#define D6 0.1144

#define A1 0
#define A2 0.329
#define A3 0.3115
#define A4 0
#define A5 0
#define A6 0


#define ALPHA1 -90
#define ALPHA2 0
#define ALPHA3 0 
#define ALPHA4 90
#define ALPHA5 90
#define ALPHA6 0

#define PI 3.141592654
#define DEG2RAD 0.01745329252
#define RAD2DEG 57.29577951

using namespace std;

void help()
{
    printf("\n");
	printf("* PLease enter the argument in the following format : \n");
	printf("* option : IJ = inverse jacobian\n");
	printf("*          IK = inverse kinematics\n");
	printf("*          FJ = forward jacobian\n");
	printf("*          FK = forward kinematics\n");
	printf("*          FJG = forward jacobian with gripper\n");
	printf("*          FKG = forward kinematics with gripper\n");
	printf("* IJ: ./tm_jacobian IJ q1 q2 q3 q4 q5 q6 x' y' z' a' b' c' \n");
	printf("* IK: ./tm_jacobian IK x y z a b c\n");
	printf("* FJ: ./tm_jacobian FJ q1 q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6 \n");
	printf("* FK: ./tm_jacobian FK q1 q2 q3 q4 q5 q6 \n");
	printf("* Input  dim: joint(radius), cartesian(m)\n");
	printf("* Output dim: joint(radius), cartesian(m)\n");
	printf("\n");
}


int main(int argc, char const *argv[])
{
	double *q_inv, *q_forward, *T;

	q_inv = new double[48];
	q_forward = new double[6];
	T = new double[16]; 

	Eigen::Matrix<float, 6, 1> jointspd = Eigen::Matrix<float, 6, 1>::Zero();
	Eigen::Matrix<float, 6, 1> effspd = Eigen::Matrix<float, 6, 1>::Zero();
	Eigen::Matrix<float, 6, 1> q_AfterHomeOffset, q_BeforeHomeOffset;
	Eigen::Matrix<float, 6, 1> home;


	home << 0, -PI*0.5, 0, PI*0.5, 0, 0;
/*
	double *q_test = new double [2];
	q_test[0] = 45*DEG2RAD;
	q_test[1] = 45*DEG2RAD;

	tm_jacobian::Forward_Kinematics_3(q_test, T);
	cout << ">>>> forward T" << endl;
	tm_jacobian::printMatrix(T,4,16);

	delete [] q_test;	
*/

	if(strncmp(argv[1],"--help",6) == 0)
	{
		help();
		return 0;
	}


	q_BeforeHomeOffset << 	strtod(argv[2],NULL),
							strtod(argv[3],NULL),
							strtod(argv[4],NULL),
							strtod(argv[5],NULL),
							strtod(argv[6],NULL),
							strtod(argv[7],NULL);

	q_AfterHomeOffset = q_BeforeHomeOffset + home;
	

	if(strncmp(argv[1],"IJ",2) == 0) // inverse jacobian 
	{
		effspd << 			strtod(argv[8],NULL),
							strtod(argv[9],NULL),
							strtod(argv[10],NULL),
							strtod(argv[11],NULL),
							strtod(argv[12],NULL),
							strtod(argv[13],NULL);

		Eigen::Matrix<float, 6, 6> Inverse_Jacobian = tm_jacobian::Inverse_Jacobian(q_AfterHomeOffset);
		jointspd = Inverse_Jacobian*effspd;
		cout << ">>>> Inverse jacobian" << endl;
		tm_jacobian::printMatrix(Inverse_Jacobian);

		cout << ">>>> joint speed" << endl;
		tm_jacobian::printMatrix(jointspd);
	}
	else if(strncmp(argv[1],"IK",2) == 0)    // inverse kinematics 
	{
		Eigen::Matrix<float, 6, 1>CartesianPosition = q_BeforeHomeOffset;
		Eigen::Matrix<float, 4, 4> T_;
	    Eigen::AngleAxisf rollAngle (CartesianPosition(3), Eigen::Vector3f::UnitZ());
	    Eigen::AngleAxisf yawAngle  (CartesianPosition(4), Eigen::Vector3f::UnitY());
	    Eigen::AngleAxisf pitchAngle(CartesianPosition(5), Eigen::Vector3f::UnitX());
	    Eigen::Quaternion<float> quaternion_matrix = rollAngle * yawAngle * pitchAngle;
	    Eigen::Matrix<float,3,3> RotationMatrix    = quaternion_matrix.matrix();
	    int num_sol;
	    
	    T_ <<   0., 0., 0., CartesianPosition(0),
	            0., 0., 0., CartesianPosition(1),
	            0., 0., 0., CartesianPosition(2),
	            0., 0., 0., 1.;
	    T_.block<3,3>(0,0) = RotationMatrix.block<3,3>(0,0);

	    if(strncmp(argv[1],"IKG",3) == 0)
	    {
		    Eigen::Matrix<float,4,4> T67, T06;
			T67 <<  1, 0, 0, 0,
					0, 1, 0, 0,
					0, 0, 1, 0.235,
					0, 0, 0, 1;

			T06 = T_*T67.inverse();
		    tm_jacobian::Matrix2DoubleArray(T06,T);
		   	num_sol =  tm_kinematics::inverse(T, q_inv);
		   	cout << ">>>> T07 " << endl;
	    	tm_jacobian::printMatrix(T_);
		}
		else
		{
			tm_jacobian::Matrix2DoubleArray(T_,T);
			num_sol =  tm_kinematics::inverse(T, q_inv);
		}
	   	
	    cout << ">>>> T06 " << endl;
	    tm_jacobian::printMatrix(T,4,16);
	    cout << ">>>> inverse q number of sols T06 : " << num_sol << endl;
		tm_jacobian::printMatrix(q_inv, 6, 6*(num_sol));
	}
	else if(strncmp(argv[1],"FJ",2) == 0)    // forward jacobian 
	{
		if(strncmp(argv[1],"FJG",3) == 0)
		{
			Eigen::Matrix<double, 6, 1> JointVelocity    = Eigen::Matrix<double, 6, 1>::Zero();
			Eigen::Matrix<double, 6, 1> EFFVelocity      = Eigen::Matrix<double, 6, 1>::Zero();
			Eigen::Matrix<double, 6, 6> Jacobian_gripper = tm_jacobian::Forward_Jacobian_gripper(q_AfterHomeOffset);
			JointVelocity << 	strtod(argv[8],NULL),
								strtod(argv[9],NULL),
								strtod(argv[10],NULL),
								strtod(argv[11],NULL),
								strtod(argv[12],NULL),
								strtod(argv[13],NULL);
			EFFVelocity = Jacobian_gripper*JointVelocity;

			cout << ">>>> jacobian" << endl;
			tm_jacobian::printMatrixd(Jacobian_gripper);
			cout << ">>>> effspd speed" << endl;
			tm_jacobian::printMatrixd(EFFVelocity);
		}
		else
		{
			jointspd << 		strtod(argv[8],NULL),
								strtod(argv[9],NULL),
								strtod(argv[10],NULL),
								strtod(argv[11],NULL),
								strtod(argv[12],NULL),
								strtod(argv[13],NULL);

			Eigen::Matrix<float, 6, 6> Jacobian = tm_jacobian::Forward_Jacobian(q_AfterHomeOffset);
			effspd = Jacobian*jointspd;
			cout << ">>>> jacobian" << endl;
			tm_jacobian::printMatrix(Jacobian);

			cout << ">>>> effspd speed" << endl;
			tm_jacobian::printMatrix(effspd);
		}
	}
	else if(strncmp(argv[1],"FK",2) == 0)    // forward kinematics 
	{
		if(strncmp(argv[1],"FKG",3) == 0)
		{
			Eigen::Matrix<double,4,4> T07, T67, T06;
			T67 <<  1, 0, 0, 0,
					0, 1, 0, 0,
					0, 0, 1, 0.235,
					0, 0, 0, 1;

			tm_jacobian::Matrix2DoubleArray(q_BeforeHomeOffset, q_forward);
			tm_jacobian::Forward_Kinematics_gripper(q_forward,T);	

			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					T07(i,j) = T[4*i+j];

			T06 = T07*T67.inverse();

			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					T[4*i+j] = T06(i,j);

			int num_sol = tm_kinematics::inverse(T, q_inv, q_forward);

			cout << ">>>> T07 " << endl;
	    	tm_jacobian::printMatrixd(T07);
	    	cout << ">>>> T06 " << endl;
	    	tm_jacobian::printMatrixd(T06);
			cout << ">>>> inverse q number of sols for T06 : " << num_sol << endl;
			tm_jacobian::printMatrix(q_inv, 6, 6*(num_sol));
		}
		else
		{
			tm_jacobian::Matrix2DoubleArray(q_BeforeHomeOffset, q_forward);
			tm_kinematics::forward(q_forward, T);
			int num_sol = tm_kinematics::inverse(T, q_inv, q_forward);

			cout << ">>>> T06" << endl;
			tm_jacobian::printMatrix(T,4,16);	
			cout << ">>>> inverse q number of sols : " << num_sol << endl;
			tm_jacobian::printMatrix(q_inv, 6, 6*(num_sol));
		}
	}
	else if(strncmp(argv[1],"TT",2) == 0)
	{
		double *q_test = new double [3];
		q_test[0] = 0;
		q_test[1] = 0;
		q_test[2] = 0;

		tm_jacobian::Forward_Kinematics_3(q_test, T);
		cout << ">>>> forward T" << endl;
		tm_jacobian::printMatrix(T,4,16);
		
		delete [] q_test;

	}

	delete [] T;
	delete [] q_inv;
	delete [] q_forward;
	return 0;
}
