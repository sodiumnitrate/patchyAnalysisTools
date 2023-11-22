#include "vec3.hpp"
#include <gtest/gtest.h>
#include <math.h>

TEST(Vec3Tests, test_assignment){
    Vec3 vec(3, 2, 5);
    ASSERT_DOUBLE_EQ(3, vec.get_x());
    ASSERT_DOUBLE_EQ(2, vec.get_y());
    ASSERT_DOUBLE_EQ(5, vec.get_z());
}

TEST(Vec3Tests, test_norm){
    Vec3 vec(2,2,5);
    ASSERT_FALSE(vec.is_normal());

    Vec3 vec2(0,0,0);
    ASSERT_FALSE(vec2.is_normal());

    Vec3 vec3(1, 0, 0);
    ASSERT_TRUE(vec3.is_normal());

    Vec3 vec4(0, -1, 0);
    ASSERT_TRUE(vec4.is_normal());
}

TEST(Vec3Tests, normalize_nonzero){
    Vec3 vec(2,2,5);
    vec.normalize();
    ASSERT_TRUE(vec.is_normal());

    Vec3 vec2(0,0,1);
    vec2.normalize();
    ASSERT_TRUE(vec2.is_normal());

    Vec3 vec3(0,0,-1);
    vec3.normalize();
    ASSERT_TRUE(vec3.is_normal());
}

TEST(Vec3Tests, normalize_zero){
    Vec3 vec(0,0,0);
    vec.normalize();
    ASSERT_FALSE(vec.is_normal());
}

TEST(Vec3Tests, get_norm_sq){
    Vec3 vec(1,2,3);
    ASSERT_DOUBLE_EQ(vec.get_norm_sq(), 14);

    Vec3 vec2(0,0,0);
    ASSERT_DOUBLE_EQ(vec2.get_norm_sq(), 0);

    Vec3 vec3(-1,0,0);
    ASSERT_DOUBLE_EQ(vec3.get_norm_sq(), 1);
}

TEST(Vec3Tests, diff){
    Vec3 vec(1,2,3);
    Vec3 vec2(1,1,1);
    Vec3 res;
    res = vec.diff(vec2);

    ASSERT_DOUBLE_EQ(res.get_x(), 0);
    ASSERT_DOUBLE_EQ(res.get_y(), 1);
    ASSERT_DOUBLE_EQ(res.get_z(), 2);
}

TEST(Vec3Tests, get_flipped){
    Vec3 vec(-1, 2, -3);
    Vec3 flipped;
    flipped = vec.get_flipped();

    ASSERT_DOUBLE_EQ(flipped.get_x(), 1);
    ASSERT_DOUBLE_EQ(flipped.get_y(), -2);
    ASSERT_DOUBLE_EQ(flipped.get_z(), 3);
}

TEST(Vec3Tests, dot){
    Vec3 vec1(1, 1, 0);
    Vec3 vec2(0, 0, 1);
    Vec3 vec3(-1, -1, 0);
    
    ASSERT_DOUBLE_EQ(vec1.dot(vec2), 0);
    ASSERT_DOUBLE_EQ(vec1.dot(vec3), -2);
}

TEST(Vec3Tests, rotate_by_euler){
    Vec3 vec(1,1,1);
    double angles[3] = {M_PI, 0, 0};
    vec.rotate(angles);
    ASSERT_DOUBLE_EQ(vec.get_x(), -1);
    ASSERT_DOUBLE_EQ(vec.get_y(), -1);
    ASSERT_DOUBLE_EQ(vec.get_z(), 1);
}

TEST(Vec3Tests, apply_pbc){
    Vec3 cell(2,2,2);
    Vec3 pos1(0.5, 0, 0);
    Vec3 pos2(2, 0, 0);

    Vec3 disp = pos2.diff(pos1);
    disp.apply_pbc(cell);

    ASSERT_DOUBLE_EQ(disp.get_x(), -0.5);
    ASSERT_DOUBLE_EQ(disp.get_y(), 0);
    ASSERT_DOUBLE_EQ(disp.get_z(), 0);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}