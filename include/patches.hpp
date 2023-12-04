/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#pragma once
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Patches{
    int n_patch = 0;
    std::vector<double> eps_vals;
    std::vector<double> lambda_vals;
    std::vector<double> cos_delta_vals;
    std::vector<Eigen::Vector3d> patch_vectors;
    std::vector<int> types;

    std::vector<std::vector<bool> > adjacency;
    double max_lambda = 0;
public:
    // constructor
    Patches();

    // add patch
    void add_patch(double eps, double lambda, double cos_delta, Eigen::Vector3d vec, int type);
    void add_patch(double eps, double lambda, double cos_delta, std::vector<double> vec, int type);

    // read from file
    void read_from_file(std::string file_name);

    // adjacency matrix
    void make_all_patches_adjacent();
    void turn_off_all_interactions();
    void make_ij_interact(int i, int j);
    bool are_ij_adjacent(int i, int j);

    // getters
    int get_n_patch();
    double get_max_lambda();
    int get_type(int i);
    Eigen::Vector3d get_vector(int i);
    double get_cos_delta(int i);
};