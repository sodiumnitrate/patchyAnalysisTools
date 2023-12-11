#include "grid.hpp"
#include <gtest/gtest.h>

TEST(GridTests, init){
    Grid<int> blah(3,3,3,0);

    Grid<float> blah_2(3,3,3,5.2);

    Grid<double> blah_3(3,3,5,23.1);
}

TEST(GridTests, reduce){
    Grid<int> blah(3,3,3,1);

    std::vector<std::vector<int> > plane = blah.sum(0);

    ASSERT_EQ(plane.size(), 3);
    for(int i = 0; i < 3; i++) ASSERT_EQ(plane[i].size(), 3);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            ASSERT_EQ(plane[i][j], 3);
        }
    }
}

TEST(GridTests, reduce_2){
    Grid<double> blah(3,5,8,0.5);

    std::vector<std::vector<double> > plane = blah.sum(2);

    ASSERT_EQ(plane.size(), 3);
    for(int i = 0; i < 3; i++) ASSERT_EQ(plane[i].size(), 5);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 5; j++){
            ASSERT_DOUBLE_EQ(plane[i][j], 0.5*8);
        }
    }

    std::vector<std::vector<double> > plane_2 = blah.sum(1);
    ASSERT_EQ(plane_2.size(), 3);
    for(int i = 0; i < 3; i++) ASSERT_EQ(plane_2[i].size(), 8);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 8; j++){
            ASSERT_DOUBLE_EQ(plane_2[i][j], 0.5*5);
        }
    }
}

TEST(GridTests, getters){
    Grid<double> blah(3,5,8,0.5);

    ASSERT_DOUBLE_EQ(blah(0,0,0), 0.5);
    ASSERT_DOUBLE_EQ(blah(2,3,7), 0.5);

    std::vector<int> size = blah.get_size();
    ASSERT_EQ(size[0], 3);
    ASSERT_EQ(size[1], 5);
    ASSERT_EQ(size[2], 8);
}

TEST(GridTests, add_subgrid){
    Grid<int> blah(3,5,8,1);
    Grid<int> subgrid(2,2,2,1);

    blah.add_subgrid(subgrid, 1, 1, 1);
    ASSERT_EQ(blah(1,1,1), 2);
    ASSERT_EQ(blah(1,1,2), 2);
    ASSERT_EQ(blah(0,0,0), 1);
    ASSERT_EQ(blah(2,2,2), 2);
}

TEST(GridTests, add_subgrid_2){
    Grid<int> blah(3,5,8,1);
    Grid<int> subgrid(2,2,2,1);

    blah.add_subgrid(subgrid, 1, 4, 1);
    ASSERT_EQ(blah(1,4,1), 2);
    ASSERT_EQ(blah(1,3,1), 1);
    ASSERT_EQ(blah(1,0,1), 2);
    ASSERT_EQ(blah(1,1,1), 1);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}