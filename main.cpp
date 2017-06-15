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
#include <Eigen/Dense>
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

	q_BeforeHomeOffset << 	strtod(argv[1],NULL),
							strtod(argv[2],NULL),
							strtod(argv[3],NULL),
							strtod(argv[4],NULL),
							strtod(argv[5],NULL),
							strtod(argv[6],NULL);
	q_BeforeHomeOffset *= DEG2RAD;
	q_AfterHomeOffset = q_BeforeHomeOffset + home;

	Eigen::Matrix<float, 6, 6> jacobian = tm_jacobian::Forward_Jacobian(q_AfterHomeOffset);
	cout << ">>>> jacobian" << endl;
	tm_jacobian::printMatrix(jacobian);


	tm_jacobian::Matrix2DoubleArray(q_BeforeHomeOffset, q_forward);
	cout << ">>>> Input q : " << endl;
	for (int i = 0; i < 6; ++i)
	{
		printf("%10.4lf ",q_forward[i]*RAD2DEG );
	}
	printf("\n");

	tm_kinematics::forward(q_forward, T);

	int num_sol = tm_kinematics::inverse(T, q_inv, q_forward);


	cout << ">>>> forward T" << endl;
	tm_jacobian::printMatrix(T,4,16);	

	cout << ">>>> inverse q number of sols : " << num_sol << endl;
	tm_jacobian::printMatrix(q_inv, 6, 6*(num_sol));

	return 0;
}
