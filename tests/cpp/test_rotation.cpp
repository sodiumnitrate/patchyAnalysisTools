#include "rotation.hpp"
#include <gtest/gtest.h>
#include <math.h>
#include <iostream>

TEST(RotationTests, init){
    Rotation rot1(0, 0, 0);

    std::vector<std::vector<double> > R = rot1.get_rotation_matrix();
    ASSERT_DOUBLE_EQ(R[0][0], 1);
    ASSERT_DOUBLE_EQ(R[1][1], 1);
    ASSERT_DOUBLE_EQ(R[2][2], 1);
    ASSERT_DOUBLE_EQ(R[0][1], 0);
    ASSERT_DOUBLE_EQ(R[0][2], 0);
    ASSERT_DOUBLE_EQ(R[1][0], 0);
    ASSERT_DOUBLE_EQ(R[1][2], 0);
    ASSERT_DOUBLE_EQ(R[2][0], 0);
    ASSERT_DOUBLE_EQ(R[2][1], 0);
}

TEST(RotationTests, identity_quat){
    Rotation rot1(0, 0, 0);
    std::vector<double> quat = rot1.get_quaternion();
    ASSERT_DOUBLE_EQ(quat[0], 0);
    ASSERT_DOUBLE_EQ(quat[1], 0);
    ASSERT_DOUBLE_EQ(quat[2], 0);
    ASSERT_DOUBLE_EQ(quat[3], 1); // real part
}

TEST(RotationTests, quaternion2){

    Rotation rot2(1, 5, 20, 3);
    std::vector<double> quat = rot2.get_quaternion();

    std::vector<std::vector<double> > R = rot2.get_rotation_matrix();
    std::vector<double> angles = rot2.get_zxz_ccw_angles();

    ASSERT_NEAR(R[0][0], -0.95402299, 0.00001);
    ASSERT_NEAR(R[0][1], -0.25287356, 0.00001);
    ASSERT_NEAR(R[0][2], 0.16091954, 0.00001);
    ASSERT_NEAR(R[1][0], 0.29885057, 0.00001);
    ASSERT_NEAR(R[1][1], -0.84367816, 0.00001);
    ASSERT_NEAR(R[1][2], 0.44597701, 0.00001);
    ASSERT_NEAR(R[2][0], 0.02298851, 0.00001);
    ASSERT_NEAR(R[2][1], 0.47356322, 0.00001);
    ASSERT_NEAR(R[2][2], 0.88045977, 0.00001);

    Rotation rot3(angles[0], angles[1], angles[2]);
    std::vector<std::vector<double> > R2 = rot3.get_rotation_matrix();
    for(int i = 0; i<3; i++){
        for(int j=0; j<3; j++){
            ASSERT_DOUBLE_EQ(R[i][j], R2[i][j]);
        }
    }

    std::vector<double> quat2 = rot3.get_quaternion();
    for(int i = 0; i < 4; i++){
        ASSERT_NEAR(quat[i], quat2[i], 0.00001);
    }
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}