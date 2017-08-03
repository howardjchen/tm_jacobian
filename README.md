# tm_jacobian : kinematics package for tm5
This is a new kinematics package for tm5 from Techman Robot.


```
* PLease enter the argument in the following format :
* option : IJ = inverse jacobian
*          IK = inverse kinematics
*          FJ = forward jacobian
*          FK = forward kinematics
* IJ: ./tm_jacobian IJ q1 q2 q3 q4 q5 q6 x' y' z' a' b' c'
* IK: ./tm_jacobian IK x y z a b c
* FJ: ./tm_jacobian FJ q1 q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6
* FK: ./tm_jacobian FK q1 q2 q3 q4 q5 q6
* Input  dim: joint(radius), cartesian(m)
* Output dim: joint(radius), cartesian(m)
```