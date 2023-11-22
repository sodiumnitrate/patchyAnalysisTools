#include "patches.hpp"

// constructor
Patches::Patches(){ n_patch = 0; }

void Patches::add_patch(double eps, double lambda, double cos_delta, Vec3 vec, int type){
    eps_vals.push_back(eps);
    lambda_vals.push_back(lambda);
    cos_delta_vals.push_back(cos_delta);
    types.push_back(type);
    vec.normalize();
    patch_vectors.push_back(vec);
    n_patch += 1;
    if (lambda > max_lambda) max_lambda = lambda;
}

void Patches::add_patch(double eps, double lambda, double cos_delta, std::vector<double> vec, int type){
    eps_vals.push_back(eps);
    lambda_vals.push_back(lambda);
    cos_delta_vals.push_back(cos_delta);
    types.push_back(type);
    Vec3 vec_(vec[0], vec[1], vec[2]);
    vec_.normalize();
    patch_vectors.push_back(vec_);
    n_patch += 1;
}

void Patches::make_all_patches_adjacent(){
    adjacency.clear();
    for(int i = 0; i < n_patch; i ++){
        std::vector<bool> row;
        for(int j = 0; j <n_patch; j++){
            row.push_back(true);
        }
        adjacency.push_back(row);
    }
}

void Patches::turn_off_all_interactions(){
    adjacency.clear();
    for(int i = 0; i < n_patch; i ++){
        std::vector<bool> row;
        for(int j = 0; j <n_patch; j++){
            row.push_back(false);
        }
        adjacency.push_back(row);
    }
}

void Patches::set_interactions(std::vector<std::vector<int> > interactions){
    this->turn_off_all_interactions();
    for(int i = 0; i < interactions.size(); i++){
        for(auto t : interactions[i]){
            adjacency[i][t] = true;
        }
    }
}

void Patches::make_ij_interact(int i, int j){
    adjacency[i][j] = true;
    adjacency[j][i] = true;
}

void Patches::read_from_file(std::string file_name){
    std::ifstream in_file(file_name);
    std::string line, ns, dummy, pm, xs, ys, zs;

    double eps, lambda, cos_delta, x, y ,z;
    int type;
    std::vector<int> interacts_with;

    if(!in_file.is_open()){
        throw "failed to open file";
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
        Vec3 vec(x,y,z);
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
    if (ct != n_patch) throw "the number of patches read does not match the declared number.";
}

int Patches::get_n_patch(){return n_patch;}

double Patches::get_max_lambda(){return max_lambda;}

bool Patches::check_type(int type, int patch_idx){
    if(types[patch_idx] == type) return true;
    return false;
}

bool Patches::are_ij_adjacent(int i, int j){
    if (adjacency[i][j]) return true;
    return false;
}

Vec3 Patches::get_vector(int i){
    return patch_vectors[i];
}

double Patches::get_cos_delta(int i){
    return cos_delta_vals[i];
}