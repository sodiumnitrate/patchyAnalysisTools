#include "rotation.hpp"

// TODO: make quaternions default, avoid extra calculations

// conversions
void Rotation::rotation_matrix_to_angles(){
    // theta
    zxz_ccw_angles[1] = std::acos(rotation_matrix[2][2]);

    if ( abs(abs(rotation_matrix[2][2]) - 1) < 1e-10 ){
        // cos(theta) = +/- 1, sin(theta) = 0
        // take phi = 0
        zxz_ccw_angles[0] = 0;
        // psi
        zxz_ccw_angles[2] = std::atan2(rotation_matrix[1][0], rotation_matrix[0][0]);
    }
    else{
        // phi
        zxz_ccw_angles[0] = std::atan2(rotation_matrix[0][2], -1*rotation_matrix[1][2]);
        // psi
        zxz_ccw_angles[2] = std::atan2(rotation_matrix[2][0], rotation_matrix[2][1]);
    }
}
void Rotation::rotation_matrix_to_quaternion(){
    double qw, qx, qy, qz, S;

    if ( rotation_matrix[2][2] < 0){
        if (rotation_matrix[0][0] > rotation_matrix[1][1]){
            S = 1 + rotation_matrix[0][0] - rotation_matrix[1][1] - rotation_matrix[2][2];
            qx = S;
            qy = rotation_matrix[0][1] + rotation_matrix[1][0];
            qz = rotation_matrix[2][0] + rotation_matrix[0][2];
            qw = rotation_matrix[1][2] - rotation_matrix[2][1];
        }
        else{
            S = 1 - rotation_matrix[0][0] + rotation_matrix[1][1] - rotation_matrix[2][2];
            qx = rotation_matrix[0][1] + rotation_matrix[1][0];
            qy = S;
            qz = rotation_matrix[1][2] + rotation_matrix[2][1];
            qw = rotation_matrix[2][0] - rotation_matrix[0][2];
        }
    }
    else{
        if ( rotation_matrix[0][0] < -1*rotation_matrix[1][1]){
            S = 1 - rotation_matrix[0][0] - rotation_matrix[1][1] + rotation_matrix[2][2];
            qx = rotation_matrix[2][0] + rotation_matrix[0][2];
            qy = rotation_matrix[1][2] + rotation_matrix[2][1];
            qz = S;
            qw = rotation_matrix[0][1] - rotation_matrix[1][0];
        }
        else{
            S = 1 + rotation_matrix[0][0] + rotation_matrix[1][1] + rotation_matrix[2][2];
            qx = rotation_matrix[1][2] - rotation_matrix[2][1];
            qy = rotation_matrix[2][0] - rotation_matrix[0][2];
            qz = rotation_matrix[0][1] - rotation_matrix[1][0];
            qw = S;
        }
    }
    
    quaternion[3] = qw * 0.5 / std::sqrt(S);
    quaternion[0] = qx * 0.5 / std::sqrt(S);
    quaternion[1] = qy * 0.5 / std::sqrt(S);
    quaternion[2] = qz * 0.5 / std::sqrt(S);
}

void Rotation::quaternion_to_rotation_matrix(){
    // qx, qy, qz, qw
    double q0, q1, q2, q3;
    q1 = quaternion[0];
    q2 = quaternion[1];
    q3 = quaternion[2];
    q0 = quaternion[3];

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
    rotation_matrix[0][1] = -1*std::sin(phi) * std::cos(theta) * std::cos(psi) - std::cos(phi) * std::sin(psi);
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
    double norm = std::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    quaternion[0] = q0;
    quaternion[1] = q1;
    quaternion[2] = q2;
    quaternion[3] = q3;

    if(abs(norm) > 1e-16){
        for(int i = 0; i < 4; i++){
            quaternion[i] /= norm;
        }
    }


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

std::vector<std::vector<double> > Rotation::get_rotation_matrix(){
    std::vector<std::vector<double> > res;
    for(int i = 0; i < 3; i++){
        std::vector<double> row;
        for (int j = 0; j < 3; j++){
            row.push_back(rotation_matrix[i][j]);
        }
        res.push_back(row);
    }
    return res;
}