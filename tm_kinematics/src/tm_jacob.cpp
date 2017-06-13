/*********************************************************************
 *                      Apache License
 *                 Version 2.0, January 2004
 *               http://www.apache.org/licenses/
 *
 * tm_jacob.cpp
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
/*
 *  tm_jacob.cpp
 *
 *  Created on: June 13, 2017
 *      Author: Howard Chen
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>

#include <Eigen/Dense>

#include "tm_kinematics/tm_jacob.h"

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

namespace tm_jacobian {

	Eigen::Matrix<float, 6, 6> Forward_Jacobian(Eigen::Matrix<float, 6,1> q)
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

	Eigen::Matrix<float, 6, 6> Inverse_Jacobian(Eigen::Matrix<float, 6,1> q)
	{
		Eigen::Matrix<float, 6, 6> jacobian = Eigen::Matrix<float, 6, 6>::Zero();
		jacobian << D6*(cos(q(0))*cos(q(4)) + sin(q(4))*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))))) - D4*cos(q(0)) - D5*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))) - A2*cos(q(1))*sin(q(0)) - A3*cos(q(1))*cos(q(2))*sin(q(0)) + A3*sin(q(0))*sin(q(1))*sin(q(2)),   D5*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))) - A2*cos(q(0))*sin(q(1)) - D6*sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))) - A3*cos(q(0))*cos(q(1))*sin(q(2)) - A3*cos(q(0))*cos(q(2))*sin(q(1)),   D5*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))) - A3*cos(q(0))*cos(q(1))*sin(q(2)) - A3*cos(q(0))*cos(q(2))*sin(q(1)),   D5*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))), -D6*(sin(q(0))*sin(q(4)) - cos(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))))),                                                                                                                                       				        						   0.,
	 				D5*(cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))) - D4*sin(q(0)) + D6*(cos(q(4))*sin(q(0)) + sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))))) + A2*cos(q(0))*cos(q(1)) + A3*cos(q(0))*cos(q(1))*cos(q(2)) - A3*cos(q(0))*sin(q(1))*sin(q(2)), - D5*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))) - A2*sin(q(0))*sin(q(1)) - A3*cos(q(1))*sin(q(0))*sin(q(2)) - A3*cos(q(2))*sin(q(0))*sin(q(1)), - D5*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))) - A3*cos(q(1))*sin(q(0))*sin(q(2)) - A3*cos(q(2))*sin(q(0))*sin(q(1)), - D5*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0)))),  D6*(cos(q(0))*sin(q(4)) - cos(q(4))*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))))),                                                                                                                                                      								   0.,
																					                                                                                                                                                                                                                                                                                                                                                                                       				 0.,                                                                                             					 A3*sin(q(1))*sin(q(2)) - D5*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - A3*cos(q(1))*cos(q(2)) - A2*cos(q(1)) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1)))),                                                                                   					   A3*sin(q(1))*sin(q(2)) - A3*cos(q(1))*cos(q(2)) - D5*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1)))),                                                                 				- D5*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - D6*sin(q(4))*(cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1)))),                                                     			 -D6*cos(q(4))*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))),                                                                                                                                                     								   0.,
																					                     	                                                                                                                                                                                                                                                                                                                                                                  			 0.,                                                                                                                                                                                                                                                                                                                                                  																		 -sin(q(0)),                                                                                                                                                                                                                                                                                                                            																	-sin(q(0)),                                                                                                                                                                                                                                                                  														 -sin(q(0)),                                   		 cos(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))),   cos(q(4))*sin(q(0)) + sin(q(4))*(cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(0))*cos(q(1))*sin(q(2)) + cos(q(0))*cos(q(2))*sin(q(1)))),
																					                                                                                                                                                                                                                                                                                                                                                                                       				 0.,                                                                                                                                                                                                                                                                                                                                                   																		  cos(q(0)),                                                                                                                                                                                                                                                                                                                             																	 cos(q(0)),                                                                                                                                                                                                                                                                   														  cos(q(0)),                                   		 cos(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1))) - sin(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))), - cos(q(0))*cos(q(4)) - sin(q(4))*(cos(q(3))*(sin(q(0))*sin(q(1))*sin(q(2)) - cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(3))*(cos(q(1))*sin(q(0))*sin(q(2)) + cos(q(2))*sin(q(0))*sin(q(1)))),
																					                                                                                                                                                                                                                                                                                                                                                                                       				 1.,                                                                                                                                                                                                                                                                                                                                                        																		 0.,                                                                                                                                                                                                                                                                                                                                  																		0.,                                                                                                                                                                                                                                                                        															 0.,                                                                   				 cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) - sin(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))),                                                    			   -sin(q(4))*(cos(q(3))*(cos(q(1))*sin(q(2)) + cos(q(2))*sin(q(1))) + sin(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))));
		return jacobian.inverse();                  
	}


	void printMatrix(Eigen::MatrixXf InputMatrix)
	{
		Eigen::MatrixXf InputTranspose = InputMatrix.transpose();
		short count = 0;
		int row = InputMatrix.rows();
		int col = InputMatrix.cols();

		for (int i = 0; i < row*col; ++i)
		{
			printf("%10.4f ", InputTranspose(i));
			if (count == row-1)
			{
				count = 0;
				printf("\n");
			}
			else
				count++;
		}
	}

	void printMatrix(double *InputMatrix, short col, int num)
	{
		short count = 0;
		if(num == 0)
			num = 1*col;
		
		for (int i = 0; i < num; ++i)
		{
			printf("%10.4f ", InputMatrix[i]);
			if (count == col-1)
			{
				count = 0;
				printf("\n");
			}
			else
				count++;
		}
		printf("\n");
	}

	void Matrix2DoubleArray(Eigen::MatrixXf InputMatrix, double *T)
	{
		Eigen::MatrixXf InputTranspose = InputMatrix.transpose();
		short row = InputMatrix.rows();
		short col = InputMatrix.cols();

		for (int i = 0; i < row*col; ++i)
			T[i] = InputTranspose(i);
	}

}