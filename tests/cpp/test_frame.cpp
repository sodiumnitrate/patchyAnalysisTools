#include "frame.hpp"
#include <gtest/gtest.h>
#include <math.h>
#include <iostream>

TEST(FrameTests, init){
    Eigen::Vector3d pos(1,1,1);
    Rotation rot(0,0,0);
    std::vector<Eigen::Vector3d> coords;
    std::vector<Rotation> orients;
    int N = 10;
    for(int i = 0; i < N; i++){
        coords.push_back(pos);
        orients.push_back(rot);
    }

    Eigen::Vector3d cell(2,2,2);

    Frame my_frame(coords, orients, cell);

    ASSERT_EQ(my_frame.get_N(), N);

    std::vector<double> cell_vec = my_frame.get_cell_as_list();
    for(int i = 0; i < 3; i++) ASSERT_EQ(cell_vec[i], cell(i));

    std::vector<int> types = my_frame.get_types();
    for( auto t : types){
        ASSERT_EQ(t, 0);
    }

    std::vector<std::vector<double> > angles = my_frame.get_orientations_as_list_of_angles();
    for( auto t : angles){
        for(int i = 0; i < 3; i++) ASSERT_EQ(t[i], 0);
    }
}

TEST(FrameTests, setters){
    Eigen::Vector3d pos(1,1,1);
    Rotation rot(0,0,0);
    std::vector<Eigen::Vector3d> coords;
    std::vector<Rotation> orients;
    int N = 10;
    for(int i = 0; i < N; i++){
        coords.push_back(pos);
        orients.push_back(rot);
    }

    Eigen::Vector3d cell(2,2,2);

    Frame my_frame(coords, orients, cell);

    ASSERT_EQ(my_frame.get_time_stamp(), 0);
    my_frame.set_time_stamp(128);
    ASSERT_EQ(my_frame.get_time_stamp(), 128);

    std::vector<int> types;
    for(int i = 0; i < N; i++){
        types.push_back(i);
    }
    my_frame.set_types(types);
    std::vector<int> types_g = my_frame.get_types();
    for(int i = 0; i < N; i++) ASSERT_EQ(types[i], types_g[i]);
}

TEST(FrameTests, test_interaction){
    std::vector<Eigen::Vector3d> coords;
    std::vector<Rotation> orients;
    Rotation rot(0,0,0);
    double sep = 1.05;
    for(int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                Eigen::Vector3d pos(i * sep, j * sep, k * sep);
                coords.push_back(pos);
                orients.push_back(rot);
            }
        }
    }

    Eigen::Vector3d cell(3 * sep, 3 * sep, 3 * sep);
    Frame frame(coords, orients, cell);

    // create the patch obj
    Patches patches;
    Eigen::Vector3d p0(1,0,0);
    patches.add_patch(1, 1.1, 0.92, p0, 0);
    Eigen::Vector3d p1(-1,0,0);
    patches.add_patch(1, 1.1, 0.92, p1, 0);
    Eigen::Vector3d p2(0,1,0);
    patches.add_patch(1, 1.1, 0.92, p2, 0);
    Eigen::Vector3d p3(0,-1,0);
    patches.add_patch(1, 1.1, 0.92, p3, 0);
    Eigen::Vector3d p4(0,0,1);
    patches.add_patch(1, 1.1, 0.92, p4, 0);
    Eigen::Vector3d p5(0,0,-1);
    patches.add_patch(1, 1.1, 0.92, p5, 0);

    ASSERT_TRUE(frame.does_ij_interact(0, 1, patches));
    ASSERT_TRUE(frame.does_ij_interact(0, 3, patches));
    ASSERT_TRUE(frame.does_ij_interact(0, 2, patches));
    ASSERT_FALSE(frame.does_ij_interact(0, 4, patches));
    ASSERT_TRUE(frame.does_ij_interact(0, 9, patches));
}

TEST(FrameTests, test_bond_list){
    std::vector<Eigen::Vector3d> coords;
    std::vector<Rotation> orients;
    Rotation rot(0,0,0);
    double sep = 1.05;
    for(int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                Eigen::Vector3d pos(i * sep, j * sep, k * sep);
                coords.push_back(pos);
                orients.push_back(rot);
            }
        }
    }

    Eigen::Vector3d cell(3 * sep, 3 * sep, 3 * sep);
    Frame frame(coords, orients, cell);

    // create the patch obj
    Patches patches;
    Eigen::Vector3d p0(1,0,0);
    patches.add_patch(1, 1.1, 0.92, p0, 0);
    Eigen::Vector3d p1(-1,0,0);
    patches.add_patch(1, 1.1, 0.92, p1, 0);
    Eigen::Vector3d p2(0,1,0);
    patches.add_patch(1, 1.1, 0.92, p2, 0);
    Eigen::Vector3d p3(0,-1,0);
    patches.add_patch(1, 1.1, 0.92, p3, 0);
    Eigen::Vector3d p4(0,0,1);
    patches.add_patch(1, 1.1, 0.92, p4, 0);
    Eigen::Vector3d p5(0,0,-1);
    patches.add_patch(1, 1.1, 0.92, p5, 0);

    // determine bond_list
    frame.determine_bond_list(patches);
    std::vector<std::vector<int> > bond_list = frame.get_bond_list();
    ASSERT_EQ(bond_list.size(), coords.size() * 3);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}