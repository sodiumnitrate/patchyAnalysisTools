#include "clusters.hpp"

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

void Clusters::clusters_from_interacting_pairs(int N){
    std::vector<std::vector<int> > neighbors;
    for (int i=0; i<N; i++){
        std::vector<int> curr;
        neighbors.push_back(curr);
    }

    for (auto t : bonds){
        neighbors[t[0]].push_back(t[1]);
        neighbors[t[1]].push_back(t[0]);
    }

    std::vector<bool> visited;
    visited.resize(N, false);
    std::queue<int> Q;
    int curr;
    for(int i = 0; i < N; i++){
        std::vector<int> curr_cluster;
        std::vector<std::vector<int> > curr_rel_bonds;
        if (visited[i]) continue;
        Q.push(i);
        curr_cluster.push_back(i);
        while(!Q.empty()){
            curr = Q.front();
            Q.pop();
            if ( visited[curr] ) continue;
            visited[curr] = true;
            curr_cluster.push_back(curr);
            for(auto t : neighbors[curr]){
                if(!visited[t]){
                    Q.push(t);
                    std::vector<int> bond;
                    bond.push_back(curr);
                    bond.push_back(t);
                    curr_rel_bonds.push_back(bond);
                }
            }
        }
        clusters.push_back(curr_cluster);
        relevant_bonds.push_back(curr_rel_bonds);
    }
}

// setting and getting private props
std::vector<std::vector<int> > Clusters::get_clusters() { return clusters; }
void Clusters::set_bonds(std::vector<std::vector<int> > bonds_){ bonds = bonds_; }
std::vector<std::vector<int> > Clusters::get_bonds(){ return bonds; }
std::vector<std::vector<std::vector<int> > > Clusters::get_relevant_bonds(){ return relevant_bonds; }