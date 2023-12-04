#include "clusters.hpp"
#include <gtest/gtest.h>

TEST(ClustersTests, init){
    std::vector<std::vector<int> > bond_list = {{0,1}, {1,2}, {2,0}, {3,4}};

    Clusters clusters(bond_list);

    std::vector<std::vector<int> > list_from_obj = clusters.get_bond_list();
    ASSERT_EQ(list_from_obj.size(), 4);

    std::vector<std::vector<int> > cluster_list = clusters.get_clusters();
    ASSERT_EQ(cluster_list.size(), 2);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}