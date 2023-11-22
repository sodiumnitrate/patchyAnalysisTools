#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "vec3.hpp"
#include "patches.hpp"
#include "clusters.hpp"

class Frame{
    int N; // number of particles
    int frame_num; // frame number
    std::vector<Vec3> coordinates; // vector of coordinates
    std::vector<Rotation> orientations; // vector of orientations
    std::vector<int> bonding; // vector containing number of bonds each particle is involved in
    std::vector<int> types; // vector containing type index for each particle
    Vec3 cell; // cell dimensions
    int time_stamp; // time stamp of the frame

    bool percolated;
    double energy;
    std::vector<int> bonds_calculated;

    Patches patches;
    Clusters cluster_info;

public:
    Frame(std::vector<Vec3> coordinates_, std::vector<Rotation> orientations_, Vec3 cell_);
    Frame(std::vector<std::vector<double> > coords, std::vector<std::vector<double> > orients, std::vector<double> cell_);
    void set_coordinates(std::vector<Vec3> coords);
    void set_coordinates(std::vector<std::vector<double> > coords);
    std::vector<Vec3> get_coordinates();
    std::vector<std::vector<double> > get_coordinates_as_list();
    std::vector<std::vector<double> > get_orientations_as_vector_of_angles();
    void set_orientations(std::vector<Rotation> orient);
    void set_orientations(std::vector<Vec3> orient);
    std::vector<Rotation> get_orientations();
    Vec3 get_cell();
    std::vector<double> get_cell_as_list();
    std::vector<int> get_types();

    int get_N();
    void set_frame_num(int num);
    int get_frame_num();


    void write_xyz(std::string file_name);
    void set_cluster_info();
    void set_patch_info(Patches patches_);
    void set_patch_info(std::string file_name);

    std::vector<std::vector<int> > get_list_of_interacting_pairs(double = 1e-6);
    bool check_calculated_bonds_against_bond_numbers();
    //std::vector<std::vector<double> > calculate_rdf(double binsize=0.1); // TODO: allow for selections
    //std::vector<std::vector<double> > calculate_sq(double g=30);
    //double get_bond_probability();
    //void find_percolating_clusters(double = 1.1);
    //bool is_system_percolated(double = 1.1);

    bool does_i_j_interact(int i, int j, double = 1e-6);

    Vec3 get_i_j_displacement(int i, int j);

    std::vector<double> get_i_j_displacement_as_vec(int i, int j);

    /*
    TODO:
    - write_frame_with_cluster_info
    - get_cluster_rg
    - write_biggest_cluster_xyz
    - get_type_fractions
    - write_frame_pdb
    - write_frame_gsd
    - get_coordination_number
    - get_density_distribution
    */
};