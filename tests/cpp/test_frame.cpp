#include "frame.hpp"
#include <gtest/gtest.h>

TEST(FrameTests, init){
    Vec3 point1(0,0,0);
    Vec3 point2(1,1.5,1);
    std::vector<Vec3> coords;
    coords.push_back(point1);
    coords.push_back(point2);

    std::vector<Rotation> orients;
    Rotation rot1(0,0,0,1);
    Rotation rot2(0,0,0,1);
    orients.push_back(rot1);
    orients.push_back(rot2);

    Vec3 cell(2,2,2);

    Frame my_frame(coords, orients, cell);

    ASSERT_EQ(my_frame.get_N(), 2);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}