#pragma once
#include <vector>
#include <cmath>

class Rotation{
    double rotation_matrix[3][3];
    double quaternion[4];  
    double zxz_ccw_angles[3];
public:
    Rotation(double rot_matrix[3][3]);
    Rotation(double quat[4]);
    Rotation(double angles[3]);
    void quaternion_to_rotation_matrix();
    void rotation_matrix_to_angles();
    void rotation_matrix_to_quaternion();

    Rotation(Rotation &t);

    std::vector<double> get_zxz_ccw_angles();
    std::vector<double> get_quaternion();

    std::vector<double> rotate_vec(std::vector<double> vec);
};