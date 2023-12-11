/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#pragma once
#include <Eigen/Dense>
#include <vector>
#include "rotation.hpp"
#include "patches.hpp"
#include "clusters.hpp"
#include "grid.hpp"
#include <string>
#include <fstream>
#include <math.h>

class Frame{
    int N;
    std::vector<Eigen::Vector3d> coordinates;
    std::vector<Rotation> orientations;
    std::vector<int> types;
    Eigen::Vector3d cell;
    int time_stamp;

    std::vector<std::vector<int> > bond_list;
    Clusters clusters;
public:
    // constructors
    Frame(std::vector<Eigen::Vector3d> coordinates_, std::vector<Rotation> orientations_, Eigen::Vector3d cell_);
    Frame(std::vector<std::vector<double> > coords, std::vector<std::vector<double> > orients, std::vector<double> cell_);

    // getters
    int get_N();
    int get_time_stamp();
    std::vector<Eigen::Vector3d> get_coordinates();
    std::vector<std::vector<double> > get_coordinates_as_list();
    std::vector<Rotation> get_orientations();
    std::vector<std::vector<double> > get_orientations_as_list_of_angles();
    Eigen::Vector3d get_cell();
    std::vector<double> get_cell_as_list();
    std::vector<int> get_types();
    std::vector<std::vector<int> > get_bond_list();
    Clusters get_clusters();

    // setters
    void set_time_stamp(int fn);
    void set_types(std::vector<int> ty);

    // write to file
    void write_xyz(std::string file_name);
    void write_xyz_file_handle(std::ofstream& out_file);

    // interactions
    Eigen::Vector3d get_ij_displacement(int i, int j);
    bool does_ij_interact(int i, int j, Patches& p);
    void determine_bond_list(Patches& p);
    void determine_clusters();
    void determine_percolation(Patches& p);
    bool is_percolated(std::vector<int> parts, Patches& p);
    bool is_percolated();
};