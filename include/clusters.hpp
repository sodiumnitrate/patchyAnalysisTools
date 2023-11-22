#pragma once
#include <vector>
#include <queue>

class Clusters{
    std::vector<std::vector<int> > clusters; // list of clusters
    std::vector<std::vector<int> > bonds; // list of bonds
    std::vector<std::vector<std::vector<int> > > relevant_bonds; // list of bonds relevant to each cluster

    std::vector<bool> percolated_clusters;
    bool percolated = false;
public:
    Clusters();
    std::vector<int> get_cluster_size();
    std::vector<std::vector<int> > get_clusters();
    void set_bonds(std::vector<std::vector<int> > bonds_);
    std::vector<std::vector<int> > get_bonds();
    std::vector<std::vector<std::vector<int> > > get_relevant_bonds();
    
    void clusters_from_interacting_pairs(int N);
    //void check_percolation(double max_lambda);
    //bool is_percolated();

    /*
    TODO:
    - get chain lenghts
    - get number of cycles
    */
};