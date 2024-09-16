#pragma once
#include <vector>
#include <cmath>

class Rotation{
    double rotation_matrix[3][3];
    double quaternion[4];  
    double zxz_ccw_angles[3];
public:
    //Rotation(double rot_matrix[3][3]);
    Rotation(double q0, double q1, double q2, double q3);
    Rotation(double phi, double theta, double psi);
    void quaternion_to_rotation_matrix();
    void rotation_matrix_to_angles();
    void rotation_matrix_to_quaternion();

    std::vector<double> get_zxz_ccw_angles();
    std::vector<double> get_quaternion();

    std::vector<double> rotate_vec(std::vector<double> vec);

    std::vector<std::vector<double> > get_rotation_matrix();
};