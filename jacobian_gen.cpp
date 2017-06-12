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

Eigen::Matrix<float,4,4> TFMatrix(float,float,float,float);
Eigen::Matrix<float, 6, 6> Jacobian_gen(Eigen::Matrix<float, 6,1>);
void PrintMatrix(Eigen::MatrixXf InputMatrix);

int main(int argc, char const *argv[])
{
	Eigen::Matrix<float, 6, 1> jointspd = Eigen::Matrix<float, 6, 1>::Zero();
	Eigen::Matrix<float, 6, 1> effspd = Eigen::Matrix<float, 6, 1>::Zero();
	Eigen::Matrix<float, 6, 1> q;
	Eigen::Matrix<float, 6, 1> home;

	home << 0, -PI*0.5, 0, PI*0.5, 0, 0;
	//jointspd << 0.3158, -0.0940, 0.1018, -0.0156, -0.3146, -0.0000;
	effspd << 0.00328786, 0.106561, 0.000812383, 0.000426278, -0.00778993, 0.00120008;

	q << 	strtod(argv[1],NULL),
			strtod(argv[2],NULL),
			strtod(argv[3],NULL),
			strtod(argv[4],NULL),
			strtod(argv[5],NULL),
			strtod(argv[6],NULL);

	//q *= DEG2RAD;
	q += home;

	Eigen::Matrix<float, 6, 6> jacobian = Jacobian_gen(q);
	cout << "======================================="<<endl;
	cout << ">>>> jacobian" << endl;
	PrintMatrix(jacobian);
	cout << "======================================="<<endl;

/*
	effspd = jacobian*jointspd;
	//effspd *= RAD2DEG;
	cout << "==================="<<endl;
	cout << ">>>> effspd speed" << endl;
	cout << effspd << endl;
	cout << "==================="<<endl;

	
	jointspd = jacobian.inverse()*effspd;
	//effspd *= RAD2DEG;
	cout << "==================="<<endl;
	cout << ">>>> jointspd speed" << endl;
	cout << jointspd << endl;
	cout << "==================="<<endl;
*/


	Eigen::Matrix<float, 4,4> T_, T_trans;
	T_  = TFMatrix(q(0),D1,ALPHA1*DEG2RAD,A1);
	T_ *= TFMatrix(q(1),D2,ALPHA2*DEG2RAD,A2);
	T_ *= TFMatrix(q(2),D3,ALPHA3*DEG2RAD,A3);
	T_ *= TFMatrix(q(3),D4,ALPHA4*DEG2RAD,A4);
	T_ *= TFMatrix(q(4),D5,ALPHA5*DEG2RAD,A5);
	T_ *= TFMatrix(q(5),D6,ALPHA6*DEG2RAD,A6);

	double InverseKinematics_T[15];
	

	T_trans = T_.transpose();

	cout << ">>>> T_" << endl;
	PrintMatrix(T_);


	cout << ">>>> InverseKinematics_T" << endl;
	int count = 0;
	for (int i = 0; i < 16; ++i)
	{
		InverseKinematics_T[i] = T_trans(i);
		printf("%10.4f ",InverseKinematics_T[i]);
		if (count == 3)
		{
			count = 0;
			printf("\n");
		}
		else
			count++;
	}

	return 0;
}


Eigen::Matrix<float, 6, 6> Jacobian_gen(Eigen::Matrix<float, 6,1> q)
{
	Eigen::Matrix<float, 6, 6> jacobian = Eigen::Matrix<float, 6, 6>::Zero();
	jacobian << D6*(cos(q(0))*cos(q(4)) + sin(q(4))*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))))) - D4*cos(q(0)) - D5*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))) - A2*cos(q(1))*sin(q(0)) - A3*cos(q(1))*cos(q(2))*sin(q(0)) + A3*sin(q(0))*sin(q(1))*sin(q(2)),   D5*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))) - A2*cos(q(0))*sin(q(1)) - D6*sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))) - A3*cos(q(0))*cos(q(1))*sin(q(2)) - A3*cos(q(0))*cos(q(2))*sin(q(1)),   D5*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))) - A3*cos(q(0))*cos(q(1))*sin(q(2)) - A3*cos(q(0))*cos(q(2))*sin(q(1)),   D5*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))), -D6*(sin(q(0))*sin(q(4)) - cos(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))))),                                                                                                                                       				        						   0.,
 				D5*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))) - D4*sin(q(0)) + D6*(cos(q(4))*sin(q(0)) + sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))))) + A2*cos(q(0))*cos(q(1)) + A3*cos(q(0))*cos(q(1))*cos(q(2)) - A3*cos(q(0))*sin(q(1))*sin(q(2)), - D5*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))) - A2*sin(q(0))*sin(q(1)) - A3*cos(q(1))*sin(q(0))*sin(q(2)) - A3*cos(q(2))*sin(q(0))*sin(q(1)), - D5*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))) - A3*cos(q(1))*sin(q(0))*sin(q(2)) - A3*cos(q(2))*sin(q(0))*sin(q(1)), - D5*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))),  D6*(cos(q(0))*sin(q(4)) - cos(q(4))*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))))),                                                                                                                                                      								   0.,
																				                                                                                                                                                                                                                                                                                                                                                                                       				 0.,                                                                                             					 A3*sin(q(1))*sin(q(2)) - D5*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - A3*cos(q(1))*cos(q(2)) - A2*cos(q(1)) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1)))),                                                                                   					   A3*sin(q(1))*sin(q(2)) - A3*cos(q(1))*cos(q(2)) - D5*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1)))),                                                                 				- D5*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1)))),                                                     			 -D6*cos(q(4))*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))),                                                                                                                                                     								   0.,
																				                     	                                                                                                                                                                                                                                                                                                                                                                  			 0.,                                                                                                                                                                                                                                                                                                                                                  																		 -sin(q(0)),                                                                                                                                                                                                                                                                                                                            																	-sin(q(0)),                                                                                                                                                                                                                                                                  														 -sin(q(0)),                                   		 cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))),   cos(q(4))*sin(q(0)) + sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))),
																				                                                                                                                                                                                                                                                                                                                                                                                       				 0.,                                                                                                                                                                                                                                                                                                                                                   																		  cos(q(0)),                                                                                                                                                                                                                                                                                                                             																	 cos(q(0)),                                                                                                                                                                                                                                                                   														  cos(q(0)),                                   		 cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))), - cos(q(0))*cos(q(4)) - sin(q(4))*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))),
																				                                                                                                                                                                                                                                                                                                                                                                                       				 1.,                                                                                                                                                                                                                                                                                                                                                        																		 0.,                                                                                                                                                                                                                                                                                                                                  																		0.,                                                                                                                                                                                                                                                                        															 0.,                                                                   				 cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))),                                                    			   -sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))));
	return jacobian;                  
}


Eigen::Matrix<float,4,4> TFMatrix(float th, float d, float Alpha, float a)
{
	Eigen::Matrix<float,4,4> T_;
	T_ << 	cos(th), -cos(Alpha)*sin(th),  sin(Alpha)*sin(th), a*cos(th),
     		sin(th),  cos(Alpha)*cos(th), -sin(Alpha)*cos(th), a*sin(th),
       			 0.,          sin(Alpha),          cos(Alpha),         d,
       		 	 0.,                  0.,                  0.,         1.;
    return T_;
}

void PrintMatrix(Eigen::MatrixXf InputMatrix)
{
	Eigen::MatrixXf Input_trans = InputMatrix.transpose();
	short count = 0;
	int row = InputMatrix.rows();
	int col = InputMatrix.cols();

	for (int i = 0; i < row*col; ++i)
	{
		printf("%10.4f ", Input_trans(i));
		if (count == row-1)
		{
			count = 0;
			printf("\n");
		}
		else
			count++;
	}
}
