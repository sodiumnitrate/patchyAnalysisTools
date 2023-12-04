/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/
#include "patches.hpp"

// constructor
Patches::Patches(){}

// add patch
void Patches::add_patch(double eps, double lambda, double cos_delta, Eigen::Vector3d vec, int type){
    eps_vals.push_back(eps);
    lambda_vals.push_back(lambda);
    if (lambda > max_lambda) max_lambda = lambda;
    cos_delta_vals.push_back(cos_delta);
    patch_vectors.push_back(vec);
    types.push_back(type);

    n_patch = eps_vals.size();
}
void Patches::add_patch(double eps, double lambda, double cos_delta, std::vector<double> vec, int type){
    Eigen::Vector3d v(vec[0], vec[1], vec[2]);
    this->add_patch(eps, lambda, cos_delta, v, type);
}

// read from file
void Patches::read_from_file(std::string file_name){
    eps_vals.clear();
    lambda_vals.clear();
    cos_delta_vals.clear();
    patch_vectors.clear();
    adjacency.clear();
    types.clear();
    max_lambda = 0;

    std::ifstream in_file(file_name);
    std::string line, ns, dummy, pm, xs, ys, zs;

    double eps, lambda, cos_delta, x, y, z;
    int type;
    std::vector<int> interacts_with();

    if(!in_file.is_open()){
        std::cout << "failed to open patch file " << file_name << std::endl;
        throw;
    }

    std::getline(in_file, line);
    std::getline(in_file, line);
    std::istringstream ss(line);
    ss >> ns >> dummy >> pm;
    n_patch = std::stoi(ns);
    if (std::stoi(pm) == 0) this->make_all_patches_adjacent();
    else this->turn_off_all_interactions();

    std::getline(in_file, line);
    std::getline(in_file, line);
    int ct = 0;
    while(std::getline(in_file, line)){
        std::getline(in_file, line);
        ss.str(line);
        ss.clear();
        ss >> dummy;
        eps = std::stod(dummy);
        ss >> dummy;
        lambda = std::stod(dummy);
        ss >> dummy;
        cos_delta = std::stod(dummy);
        std::getline(in_file, line);
        std::getline(in_file, line);
        ss.str(line);
        ss.clear();
        ss >> dummy;
        x = std::stod(dummy);
        ss >> dummy;
        y = std::stod(dummy);
        ss >> dummy;
        z = std::stod(dummy);
        ss >> dummy;
        type = std::stoi(dummy);
        Eigen::Vector3d vec(x,y,z);
        this->add_patch(eps, lambda, cos_delta, vec, type);
        std::getline(in_file, line);
        std::getline(in_file, line);
        ss.str(line);
        ss.clear();
        while(ss >> dummy){
            this->make_ij_interact(ct, std::stoi(dummy));
        }
        ct += 1;
    }
}

// adjacency
void Patches::make_all_patches_adjacent(){
    adjacency.clear();
    for(int i = 0; i < n_patch; i++){
        std::vector<bool> row;
        for(int j = 0; j < n_patch; j++){
            row.push_back(true);
        }
        adjacency.push_back(row);
    }
}
void Patches::turn_off_all_interactions(){
    adjacency.clear();
    for(int i = 0; i < n_patch; i++){
        std::vector<bool> row;
        for(int j = 0; j < n_patch; j++){
            row.push_back(false);
        }
        adjacency.push_back(row);
    }
}
void Patches::make_ij_interact(int i, int j){
    adjacency[i][j] = true;
    adjacency[j][i] = true;
}
bool Patches::are_ij_adjacent(int i, int j){
    return adjacency[i][j];
}

// getters
int Patches::get_n_patch(){ return n_patch; }
double Patches::get_max_lambda(){ return max_lambda; }
int Patches::get_type(int i){ return types[i]; }
Eigen::Vector3d Patches::get_vector(int i){ return patch_vectors[i]; }
double Patches::get_cos_delta(int i){ return cos_delta_vals[i]; }