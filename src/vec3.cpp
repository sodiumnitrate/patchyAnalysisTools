#include "vec3.hpp"


// constructors
Vec3::Vec3(){};

Vec3::Vec3(double x_, double y_, double z_){
    x = x_;
    y = y_;
    z = z_;
}

Vec3::Vec3(std::vector<double> pos){
    if (pos.size() != 3) throw "input vector size is not 3.";
    x = pos[0];
    y = pos[1];
    z = pos[2];
}

// methods to get components
double Vec3::get_x() { return x; }
double Vec3::get_y() { return y; }
double Vec3::get_z() { return z; }

std::vector<double> Vec3::get_vec() {
    std::vector<double> res;
    res.push_back(x);
    res.push_back(y);
    res.push_back(z);
    return res;
}

void Vec3::set_vec(std::vector<double> vec){
    x = vec[0];
    y = vec[1];
    z = vec[2];
}

// dot product
double Vec3::dot(Vec3 another_vec){
    std::vector<double> other = another_vec.get_vec();
    return other[0] * x + other[1] * y + other[2] * z;
}

void Vec3::rotate(Rotation& rot){
    std::vector<double> vec = this->get_vec();
    std::vector<double> res = rot.rotate_vec(vec);

    this->set_vec(res);
}

void Vec3::rotate(double angles[3]){
    Rotation rot(angles[0], angles[1], angles[2]);
    std::vector<double> vec = this->get_vec();
    std::vector<double> res = rot.rotate_vec(vec);
    this->set_vec(res);
}

bool Vec3::is_normal(double tol){
    double norm_sq = x * x + y * y + z * z;
    return std::abs(norm_sq - 1) >= tol;
}

void Vec3::normalize(){
    double norm_sq = x * x + y * y + z * z;
    double norm = std::sqrt(norm_sq);
    // if norm is zero, don't do anything (is this a reasonable tolerance?)
    if(std::abs(norm) < 1e-16) return;
    x /= norm;
    y /= norm;
    z /= norm;
}

Vec3 Vec3::get_flipped(){
    Vec3 new_vec(-1*x,-1*y,-1*z);
    return new_vec;
}

double Vec3::get_norm_sq(){
    return x * x + y * y + z * z;
}

Vec3 Vec3::diff(Vec3 another){
    double xn, yn, zn;
    xn = x - another.get_x();
    yn = y - another.get_y();
    zn = z - another.get_z();
    Vec3 new_vec(xn, yn, zn);
    return new_vec;
}

void Vec3::apply_pbc(Vec3 cell){
    std::vector<double> cellv = cell.get_vec();
    std::vector<double> vec = this->get_vec();

    for(int i = 0; i < 3; i++){
        double hbox = cellv[i] / 2;
        while(vec[i] > hbox) vec[i] -= cellv[i];
        while(vec[i] < -1*hbox) vec[i] += cellv[i];
    }

    x = cellv[0];
    y = cellv[1];
    z = cellv[2];
}