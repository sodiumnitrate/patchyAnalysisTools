/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#pragma once
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

class Clusters{
    std::vector<std::vector<int> > clusters; // list of clusters
    std::vector<std::vector<int> > bonds; // list of bonds
public:
    // constructor
    Clusters(std::vector<std::vector<int> > bond_list);

    // getters
    std::vector<std::vector<int> > get_clusters();
    std::vector<std::vector<int> > get_bond_list();
};