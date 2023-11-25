#include "patches.hpp"
#include <gtest/gtest.h>

TEST(PatchesTests, init){
    Patches patchlist;
    ASSERT_EQ(patchlist.get_n_patch(), 0);

    Vec3 patch_vec(1,0,0);

    patchlist.add_patch(1, 1.5, 0.92, patch_vec, 0);
    ASSERT_EQ(patchlist.get_n_patch(), 1);
    ASSERT_DOUBLE_EQ(patchlist.get_max_lambda(), 1.5);

    Vec3 vec = patchlist.get_vector(0);
    ASSERT_DOUBLE_EQ(vec.get_x(), patch_vec.get_x());
    ASSERT_DOUBLE_EQ(vec.get_y(), patch_vec.get_y());
    ASSERT_DOUBLE_EQ(vec.get_z(), patch_vec.get_z());

    ASSERT_DOUBLE_EQ(patchlist.get_cos_delta(0), 0.92);
    ASSERT_TRUE(patchlist.check_type(0, 0));
}

TEST(PatchesTests, read_from_file){
    Patches patches;
    // TODO: make this path more robust
    patches.read_from_file("../tests/aux_files/patches.dat");
    ASSERT_EQ(patches.get_n_patch(), 5);
    ASSERT_DOUBLE_EQ(patches.get_max_lambda(), 1.1);

    for(int i = 0; i <  patches.get_n_patch(); i++){
        ASSERT_DOUBLE_EQ(patches.get_cos_delta(i), 0.92);
        ASSERT_TRUE(patches.get_vector(i).is_normal());
        for(int j = 0; j < patches.get_n_patch(); j++){
            ASSERT_TRUE(patches.are_ij_adjacent(i,j));
        }
    }

    ASSERT_TRUE(patches.check_type(0, 0));
    ASSERT_TRUE(patches.check_type(0, 1));
    ASSERT_TRUE(patches.check_type(1, 2));
    ASSERT_TRUE(patches.check_type(1, 3));
    ASSERT_TRUE(patches.check_type(1, 4));   


}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}