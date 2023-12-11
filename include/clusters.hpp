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
    std::vector<bool> percolated; // whether or not each cluster is percolated
public:
    // constructor
    Clusters();
    Clusters(std::vector<std::vector<int> > bond_list);

    // getters
    std::vector<std::vector<int> > get_clusters();
    std::vector<std::vector<int> > get_bond_list();
    std::vector<bool> get_perc_list();

    // setters
    void set_percolation_info(std::vector<bool> perc);
};