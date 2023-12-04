#include "rotation.hpp"
#include <gtest/gtest.h>
#include <math.h>
#include <iostream>

TEST(RotationTests, init_identity){
    Rotation rot1(0, 0, 0);
    Eigen::Quaterniond q = rot1.get_quaternion();
    ASSERT_DOUBLE_EQ(q.x(), 0);
    ASSERT_DOUBLE_EQ(q.y(), 0);
    ASSERT_DOUBLE_EQ(q.z(), 0);
    ASSERT_DOUBLE_EQ(q.w(), 1);

    Eigen::Matrix3d rot_mat = rot1.get_rotation_matrix();
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            if(i == j) ASSERT_DOUBLE_EQ(rot_mat(i,j), 1);
            else ASSERT_DOUBLE_EQ(rot_mat(i,j), 0);
        }
    }

    std::vector<double> angles = rot1.get_zxz_ccw_angles();
    for(int i = 0; i < 3; i++) ASSERT_EQ(angles[i], 0);
}

TEST(RotationTests, init_by_quat){
    Eigen::Quaterniond q(1, 0, 0, 0);
    Rotation rot1(q);

    Eigen::Matrix3d rot_mat = rot1.get_rotation_matrix();
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            if(i == j) ASSERT_DOUBLE_EQ(rot_mat(i,j), 1);
            else ASSERT_DOUBLE_EQ(rot_mat(i,j), 0);
        }
    }

    std::vector<double> angles = rot1.get_zxz_ccw_angles();
    for(int i = 0; i < 3; i++) ASSERT_EQ(angles[i], 0);
}

TEST(RotationTests, angle_conversion){
    double phi, theta, psi;
    phi = 0.5;
    theta = 1.2;
    psi = 2.1;

    Rotation rot1(phi, theta, psi);
    std::vector<double> angles = rot1.get_zxz_ccw_angles();
    ASSERT_DOUBLE_EQ(angles[0], 0.5);
    ASSERT_DOUBLE_EQ(angles[1], 1.2);
    ASSERT_DOUBLE_EQ(angles[2], 2.1);

    Eigen::Matrix3d mat = rot1.get_rotation_matrix();
    Rotation rot2(mat);
    std::vector<double> angles2 = rot2.get_zxz_ccw_angles();
    for (int i = 0; i < 3; i++) ASSERT_DOUBLE_EQ(angles[i], angles2[i]);

    Eigen::Quaterniond q1 = rot1.get_quaternion();
    Eigen::Quaterniond q2 = rot2.get_quaternion();
    ASSERT_DOUBLE_EQ(q1.x(), q2.x());
    ASSERT_DOUBLE_EQ(q1.y(), q2.y());
    ASSERT_DOUBLE_EQ(q1.z(), q2.z());
    ASSERT_DOUBLE_EQ(q1.w(), q2.w());
}

TEST(RotationTests, theta_0){
    double phi, theta, psi;
    phi = 0.5;
    theta = 0;
    psi = 2.1;

    Rotation rot1(phi, theta, psi);
    Eigen::Matrix3d mat = rot1.get_rotation_matrix();
    Rotation rot2(mat);
    std::vector<double> angles = rot2.get_zxz_ccw_angles();
    ASSERT_DOUBLE_EQ(angles[1], 0);
    ASSERT_DOUBLE_EQ(angles[0], 0);
    ASSERT_DOUBLE_EQ(angles[2], phi + psi);
}

TEST(RotationTests, rotate_vec){
    Eigen::Vector3d vec(1,1,1);

    Rotation rot1(0,0,0);
    vec = rot1.rotate_vec(vec);
    for(int i = 0; i < 3; i++) ASSERT_DOUBLE_EQ(vec(i), 1);

    Rotation rot2(M_PI_2, 0, 0);
    vec = rot2.rotate_vec(vec);
    ASSERT_DOUBLE_EQ(vec(0), -1);
    ASSERT_DOUBLE_EQ(vec(1), 1);
    ASSERT_DOUBLE_EQ(vec(2), 1);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}