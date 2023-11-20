#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "vec3.hpp"

class Patches{
    int n_patch; // number of patches
    std::vector<double> eps_vals; // list that holds epsilon values (interaction strength)
    std::vector<double> lambda_vals; // list that holds lambda values
    std::vector<double> cos_delta_vals; // list that holds cos_delta values (patch widths)
    std::vector<Vec3> patch_vectors;
    std::vector<int> types; // which particle types the patch will be active on

    std::vector<std::vector<bool> > adjacency; // adjacency matrix for patch interactions
    double max_lambda = 0;
public:
    Patches();
    void add_patch(double eps, double lambda, double cos_delta, Vec3 vec, int type);
    void add_patch(double eps, double lambda, double cos_delta, std::vector<double> vec, int type);
    void make_all_patches_adjacent();
    void set_interactions(std::vector<std::vector<int> > interactions);
    void turn_off_all_interactions();
    void make_ij_interact(int i, int j);

    void read_from_file(std::string file_name);

    int get_n_patch();
    double get_max_lambda();
    bool check_type(int type, int patch_idx);
    bool are_ij_adjacent(int i, int j);

    Vec3 get_vector(int i);
    double get_cos_delta(int i);
};