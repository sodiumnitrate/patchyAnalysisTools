/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#include "rotation.hpp"

// conversions
void Rotation::rotation_matrix_to_angles(){
    // theta
    zxz_ccw_angles[1] = std::acos(rotation_matrix(2,2));

    if ( abs(abs(rotation_matrix(2,2)) - 1) < 1e-10){
        // cos(theta) = +/- 1, sin(theta) = 0
        // take phi = 0
        zxz_ccw_angles[0] = 0;
        // psi
        zxz_ccw_angles[2] = std::atan2(rotation_matrix(1,0), rotation_matrix(0,0));
    }
    else{
        // phi
        zxz_ccw_angles[0] = std::atan2(rotation_matrix(0,2), -1*rotation_matrix(1,2));
        // psi
        zxz_ccw_angles[2] = std::atan2(rotation_matrix(2,0), rotation_matrix(2,1)); 
    }
}
void Rotation::rotation_matrix_to_quaternion(){
    q = rotation_matrix;
}
void Rotation::quaternion_to_rotation_matrix(){
    rotation_matrix = q.toRotationMatrix();
}

// constructors
Rotation::Rotation(double phi, double theta, double psi){
    this->set_zxz_ccw_angles(phi, theta, psi);
}
Rotation::Rotation(Eigen::Quaterniond q_){
    this->set_quaternion(q_);
}
Rotation::Rotation(Eigen::Matrix3d rotation_matrix_){
    this->set_rotation_matrix(rotation_matrix_);
}

// setters
void Rotation::set_zxz_ccw_angles(double phi, double theta, double psi){
    zxz_ccw_angles[0] = phi;
    zxz_ccw_angles[1] = theta;
    zxz_ccw_angles[2] = psi;

    rotation_matrix(0,0) = std::cos(phi) * std::cos(psi) - std::cos(theta) * std::sin(phi) * std::sin(psi);
    rotation_matrix(0,1) = -1*std::sin(phi) * std::cos(theta) * std::cos(psi) - std::cos(phi) * std::sin(psi);
    rotation_matrix(0,2) = std::sin(phi) * std::sin(theta);
    rotation_matrix(1,0) = std::cos(psi) * std::sin(phi) + std::cos(phi) * std::cos(theta) * std::sin(psi);
    rotation_matrix(1,1) = std::cos(phi) * std::cos(theta) * std::cos(psi) - std::sin(phi) * std::sin(psi);
    rotation_matrix(1,2) = -1*std::cos(phi) * std::sin(theta);
    rotation_matrix(2,0) = std::sin(psi) * std::sin(theta);
    rotation_matrix(2,1) = std::cos(psi) * std::sin(theta);
    rotation_matrix(2,2) = std::cos(theta);

    this->rotation_matrix_to_quaternion();
}
void Rotation::set_quaternion(Eigen::Quaterniond q_){
    q = q_.normalized();

    this->quaternion_to_rotation_matrix();
    this->rotation_matrix_to_angles();
}
void Rotation::set_rotation_matrix(Eigen::Matrix3d rotation_matrix_){
    rotation_matrix = rotation_matrix_;

    this->rotation_matrix_to_angles();
    this->rotation_matrix_to_quaternion();
}

// getters
std::vector<double> Rotation::get_zxz_ccw_angles(){
    std::vector<double> res;
    for (int i = 0; i < 3; i++) res.push_back(zxz_ccw_angles[i]);
    return res;
}
Eigen::Quaterniond Rotation::get_quaternion(){ return q; }
Eigen::Matrix3d Rotation::get_rotation_matrix(){ return rotation_matrix; }

// rotate vector
Eigen::Vector3d Rotation::rotate_vec(Eigen::Vector3d vec){
    return q * vec;
}