#include "rotation.hpp"

// TODO: make quaternions default, avoid extra calculations

// conversions
void Rotation::rotation_matrix_to_angles(){
    zxz_ccw_angles[1] = std::acos(rotation_matrix[2][2]);
    zxz_ccw_angles[0] = std::atan2(-1*rotation_matrix[1][2], rotation_matrix[0][2]);
    zxz_ccw_angles[2] = std::atan2(rotation_matrix[2][0], rotation_matrix[2][1]);
}
void Rotation::rotation_matrix_to_quaternion(){
    double trace = 0;
    double qw, qx, qy, qz, S, invS;
    for(int i = 0; i < 3; i++) trace += rotation_matrix[i][i];

    if ( trace > 0){
        S = 2 * std::sqrt(trace + 1);
        invS = 1. / S;
        qw = 0.25 * S;
        qx = (rotation_matrix[2][1] - rotation_matrix[1][2]) * invS;
        qy = (rotation_matrix[0][2] - rotation_matrix[2][0]) * invS;
        qz = (rotation_matrix[1][0] - rotation_matrix[0][1]) * invS;
    }
    else if ((rotation_matrix[0][0] > rotation_matrix[1][1]) && (rotation_matrix[0][0] > rotation_matrix[2][2])){
        S = 2 * std::sqrt(1 + rotation_matrix[0][0] - rotation_matrix[1][1] - rotation_matrix[2][2]);
        invS = 1. / S;
        qw = (rotation_matrix[2][1] - rotation_matrix[1][2]) * invS;
        qx = 0.25 * S;
        qy = (rotation_matrix[0][1] + rotation_matrix[1][0]) * invS;
        qz = (rotation_matrix[0][2] + rotation_matrix[2][0]) * invS;
    }
    else if (rotation_matrix[1][1] > rotation_matrix[2][2]){
        S = std::sqrt(1.0 + rotation_matrix[1][1] - rotation_matrix[0][0] - rotation_matrix[2][2]) * 2; 
        invS = 1. / S;
        qw = (rotation_matrix[0][2] - rotation_matrix[2][0]) * invS;
        qx = (rotation_matrix[0][1] + rotation_matrix[1][0]) * invS; 
        qy = 0.25 * S;
        qz = (rotation_matrix[1][2] + rotation_matrix[2][1]) * invS; 
    }
    else { 
        S = sqrt(1.0 + rotation_matrix[2][2] - rotation_matrix[0][0] - rotation_matrix[1][1]) * 2;
        qw = (rotation_matrix[1][0] - rotation_matrix[0][1]) * invS;
        qx = (rotation_matrix[0][2] + rotation_matrix[2][0]) * invS;
        qy = (rotation_matrix[1][2] + rotation_matrix[2][1]) * invS;
        qz = 0.25 * S;
    }

    quaternion[0] = qw;
    quaternion[1] = qx;
    quaternion[2] = qy;
    quaternion[3] = qz;
}

void Rotation::quaternion_to_rotation_matrix(){
    double q0, q1, q2, q3;
    q0 = quaternion[0];
    q1 = quaternion[1];
    q2 = quaternion[2];
    q3 = quaternion[3];

    rotation_matrix[0][0] = 2 * (q0 * q0 + q1 * q1) - 1;
    rotation_matrix[0][1] = 2 * (q1 * q2 - q0 * q3);
    rotation_matrix[0][2] = 2 * (q1 * q3 + q0 * q2);
     
    rotation_matrix[1][0] = 2 * (q1 * q2 + q0 * q3);
    rotation_matrix[1][1] = 2 * (q0 * q0 + q2 * q2) - 1;
    rotation_matrix[1][2] = 2 * (q2 * q3 - q0 * q1);
     
    rotation_matrix[2][0] = 2 * (q1 * q3 - q0 * q2);
    rotation_matrix[2][1] = 2 * (q2 * q3 + q0 * q1);
    rotation_matrix[2][2] = 2 * (q0 * q0 + q3 * q3) - 1;
}

// constructors
Rotation::Rotation(double phi, double theta, double psi){
    // form rotation matrix, and quaternion, given 3 Euler angles in ZXZ CCW'
    zxz_ccw_angles[0] = phi;
    zxz_ccw_angles[1] = theta;
    zxz_ccw_angles[2] = psi;

    rotation_matrix[0][0] = std::cos(phi) * std::cos(psi) - std::cos(theta) * std::sin(phi) * std::sin(psi);
    rotation_matrix[0][1] = -1*std::sin(phi) * std::cos(theta) * std::cos(psi) - std::cos(phi) * std::cos(psi);
    rotation_matrix[0][2] = std::sin(phi) * std::sin(theta);
    rotation_matrix[1][0] = std::cos(psi) * std::sin(phi) + std::cos(phi) * std::cos(theta) * std::sin(psi);
    rotation_matrix[1][1] = std::cos(phi) * std::cos(theta) * std::cos(psi) - std::sin(phi) * std::sin(psi);
    rotation_matrix[1][2] = -1*std::cos(phi) * std::sin(theta);
    rotation_matrix[2][0] = std::sin(psi) * std::sin(theta);
    rotation_matrix[2][1] = std::cos(psi) * std::sin(theta);
    rotation_matrix[2][2] = std::cos(theta);

    this->rotation_matrix_to_quaternion();
}

Rotation::Rotation(double q0, double q1, double q2, double q3){
    quaternion[0] = q0;
    quaternion[1] = q1;
    quaternion[2] = q2;
    quaternion[3] = q3;

    this->quaternion_to_rotation_matrix();
    this->rotation_matrix_to_angles();
}

std::vector<double> Rotation::rotate_vec(std::vector<double> vec){
    std::vector<double> res;
    if (vec.size() != 3) throw "not a 3d vector";
    for(int i=0; i<3;i++){
        double val = 0;
        for(int j=0; j<3; j++){
            val += rotation_matrix[i][j] * vec[j];
        }
        res.push_back(val);
    }
    return res;
}


// access private props
std::vector<double> Rotation::get_zxz_ccw_angles(){
    std::vector<double> res;
    for(int i = 0; i < 3; i++) res.push_back(zxz_ccw_angles[i]);
    return res;
}

std::vector<double> Rotation::get_quaternion(){
    std::vector<double> res;
    for(int i = 0; i < 4; i++) res.push_back(quaternion[i]);
    return res;
}