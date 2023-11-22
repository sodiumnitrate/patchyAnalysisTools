#pragma once
#include <vector>
#include "rotation.hpp"

class Vec3{
    double x;
    double y;
    double z;
public:
    Vec3();
    Vec3(double x_, double y_, double z_);
    Vec3(std::vector<double> pos);
    std::vector<double> get_vec();
    void set_vec(std::vector<double> vec);
    double get_x();
    double get_y();
    double get_z();

    double dot(Vec3 another_vec);

    void rotate(Rotation& rot);
    void rotate(double angles[3]);

    bool is_normal(double = 1e-6);
    void normalize();

    Vec3 diff(Vec3 another);
    void apply_pbc(Vec3 cell);
    double get_norm_sq();
    Vec3 get_flipped();
};