#pragma once
#include <vector>

class Clusters{
    std::vector<std::vector<int> > clusters; // list of clusters
    std::vector<std::vector<int> > bonds; // list of bonds
    std::vector<std::vector<std::vector<int> > > relevant_bonds; // list of bonds relevant to each cluster

    std::vector<bool> percolated_clusters;
public:
    Clusters();
    std::vector<int> get_cluster_size();
    void set_clusters(std::vector<std::vector<int> > clusters_);
    std::vector<std::vector<int> > get_clusters();
    void set_bonds(std::vector<std::vector<int> > bonds_);
    std::vector<std::vector<int> > get_bonds();
    void set_relevant_bonds(std::vector<std::vector<std::vector<int> > > rel_bonds);
    std::vector<std::vector<std::vector<int> > > get_relevant_bonds();

    /*
    TODO:
    - get chain lenghts
    - get number of cycles
    */
};