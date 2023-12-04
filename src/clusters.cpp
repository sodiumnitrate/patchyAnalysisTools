/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#include "clusters.hpp"

Clusters::Clusters(std::vector<std::vector<int> > bond_list){
    bonds = bond_list;

    // determine clusters
    int i, j;
    std::unordered_map<int, std::vector<int> > neighbors;
    std::vector<int> parts;
    for (auto bond : bonds){
        i = bond[0];
        j = bond[1];
        if(neighbors.find(i) == neighbors.end()){
            std::vector<int> n = {j};
            neighbors[i] = n;
            parts.push_back(i);
        }
        else{
            neighbors[i].push_back(j);
        }
        if(neighbors.find(j) == neighbors.end()){
            std::vector<int> n = {i};
            neighbors[j] = n;
            parts.push_back(j);
        }
        else{
            neighbors[j].push_back(i);
        }
    }

    std::unordered_set<int> visited;
    std::queue<int> Q;
    int curr;
    for (auto part : parts){
        std::vector<int> curr_cluster;
        if (visited.find(part) != visited.end()) continue;
        visited.insert(part);
        curr_cluster.push_back(part);
        Q.push(part);
        while(!Q.empty()){
            curr = Q.front();
            Q.pop();
            if (visited.find(curr) != visited.end()) continue;
            curr_cluster.push_back(curr);
            for (auto n : neighbors[curr]){
                if(visited.find(n) == visited.end()){
                    Q.push(n);
                }
            }
        }
        clusters.push_back(curr_cluster);
    }
}

// getters
std::vector<std::vector<int> > Clusters::get_clusters(){ return clusters; }
std::vector<std::vector<int> > Clusters::get_bond_list(){ return bonds; }