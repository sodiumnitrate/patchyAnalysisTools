#include "include/clusters.hpp"

// constructor
Clusters::Clusters(){};

// get cluster sizes
std::vector<int> Clusters::get_cluster_size(){
    std::vector<int> res;
    for (auto t : clusters){
        res.push_back(t.size());
    }
    return res;
}

// setting and getting private props
void Clusters::set_clusters(std::vector<std::vector<int> > clusters_){ clusters = clusters_; }
std::vector<std::vector<int> > Clusters::get_clusters() { return clusters; }
void Clusters::set_bonds(std::vector<std::vector<int> > bonds_){ bonds = bonds_; }
std::vector<std::vector<int> > Clusters::get_bonds(){ return bonds; }
void Clusters::set_relevant_bonds(std::vector<std::vector<std::vector<int> > > rel_bonds){ relevant_bonds = rel_bonds; }
std::vector<std::vector<std::vector<int> > > Clusters::get_relevant_bonds(){ return relevant_bonds; }