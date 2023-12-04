#include "patches.hpp"
#include <gtest/gtest.h>
#include <math.h>
#include <iostream>

TEST(PatchesTests, init){
    Patches my_patches;
    ASSERT_EQ(my_patches.get_n_patch(), 0);
    ASSERT_DOUBLE_EQ(my_patches.get_max_lambda(), 0);
}

TEST(PatchesTests, add_patch_and_adjacency){
    Patches my_patches;
    Eigen::Vector3d vec(0,0,1);
    my_patches.add_patch(1, 1.1, 0.92, vec, 0);
    ASSERT_EQ(my_patches.get_n_patch(), 1);
    ASSERT_DOUBLE_EQ(my_patches.get_max_lambda(), 1.1);

    my_patches.add_patch(1, 1.5, 0.92, vec, 0);
    ASSERT_EQ(my_patches.get_n_patch(), 2);
    ASSERT_DOUBLE_EQ(my_patches.get_max_lambda(), 1.5);

    my_patches.add_patch(1, 1.5, 0.92, vec, 0);
    ASSERT_EQ(my_patches.get_n_patch(), 3);
    ASSERT_DOUBLE_EQ(my_patches.get_max_lambda(), 1.5);

    my_patches.turn_off_all_interactions();
    my_patches.make_ij_interact(0, 1);
    ASSERT_TRUE(my_patches.are_ij_adjacent(0, 1));
    ASSERT_FALSE(my_patches.are_ij_adjacent(0, 2));

    my_patches.make_all_patches_adjacent();
    for(int i = 0; i < my_patches.get_n_patch(); i++){
        for(int j = 0; j < my_patches.get_n_patch(); j++){
            ASSERT_TRUE(my_patches.are_ij_adjacent(i, j));
        }
    }
}

TEST(PatchesTests, getters){
    Patches my_patches;
    Eigen::Vector3d vec(0,0,1);
    my_patches.add_patch(1, 1.1, 0.92, vec, 0);
    my_patches.add_patch(1, 1.3, 0.92, vec, 0);
    my_patches.add_patch(2, 1.1, 0.92, vec, 1);

    ASSERT_DOUBLE_EQ(my_patches.get_cos_delta(0), 0.92);
    ASSERT_EQ(my_patches.get_type(0), 0);
    ASSERT_EQ(my_patches.get_type(2), 1);
    Eigen::Vector3d vec2 = my_patches.get_vector(0);
    for(int i = 0; i < 3; i++) ASSERT_DOUBLE_EQ(vec(i), vec2(i));
}

TEST(PatchesTests, read_from_file){
    Patches my_obj;
    my_obj.read_from_file("../tests/aux_files/patches.dat");
    // there should be 5 patches
    ASSERT_EQ(my_obj.get_n_patch(), 5);

    // all should be interacting with each other
    for(int i = 0; i < my_obj.get_n_patch(); i++){
        for(int j = 0; j < my_obj.get_n_patch(); j++){
            ASSERT_TRUE(my_obj.are_ij_adjacent(i, j));
        }
    }

    ASSERT_DOUBLE_EQ(my_obj.get_max_lambda(), 1.1);

    Eigen::Vector3d vec = my_obj.get_vector(0);
    ASSERT_DOUBLE_EQ(vec(0), 1);
    ASSERT_DOUBLE_EQ(vec(1), 0);
    ASSERT_DOUBLE_EQ(vec(2), 0);

    ASSERT_EQ(my_obj.get_type(0), 0);
    ASSERT_EQ(my_obj.get_type(1), 0);
    ASSERT_EQ(my_obj.get_type(2), 1);
    ASSERT_EQ(my_obj.get_type(3), 1);
    ASSERT_EQ(my_obj.get_type(4), 1);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}