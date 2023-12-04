/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#pragma once
#include <vector>
#include <cmath>
#include <Eigen/Geometry>

class Rotation{
    Eigen::Matrix3d rotation_matrix;
    Eigen::Quaterniond q;
    double zxz_ccw_angles[3];
public:
    Rotation(Eigen::Quaterniond q_);
    Rotation(double phi, double theta, double psi);
    Rotation(Eigen::Matrix3d rotation_matrix_);

    // conversions
    void quaternion_to_rotation_matrix();
    void rotation_matrix_to_angles();
    void rotation_matrix_to_quaternion();

    // getters and setters
    void set_zxz_ccw_angles(double phi_, double theta_, double psi_);
    void set_quaternion(Eigen::Quaterniond q_);
    void set_rotation_matrix(Eigen::Matrix3d rotation_matrix_);

    std::vector<double> get_zxz_ccw_angles();
    Eigen::Quaterniond get_quaternion();
    Eigen::Matrix3d get_rotation_matrix();

    // rotate vector
    Eigen::Vector3d rotate_vec(Eigen::Vector3d vec);
};