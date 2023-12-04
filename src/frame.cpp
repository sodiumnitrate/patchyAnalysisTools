/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/
#include "frame.hpp"
#include <iostream>

// constructors
Frame::Frame(std::vector<Eigen::Vector3d> coordinates_, std::vector<Rotation> orientations_, Eigen::Vector3d cell_){
    if (coordinates_.size() != orientations_.size()){
        std::cout << "ERROR: coordinates and orientations must have the same number of elements." << std::endl;
        throw;
    }
    coordinates = coordinates_;
    orientations = orientations_;
    cell = cell_;

    N = coordinates.size();

    // automatically set all types to 0
    types.resize(N);
    std::fill(types.begin(), types.end(), 0);

    // by default set frame_num to 0
    time_stamp = 0;
}
Frame::Frame(std::vector<std::vector<double> > coords, std::vector<std::vector<double> > orients, std::vector<double> cell_){
    if (coords.size() != orients.size()){
        std::cout << "ERROR: coordinates and orientations must have the same number of elements." << std::endl;
        throw;
    }

    if (cell.size() != 3){
        std::cout << "ERROR: cell must have 3 elements." << std::endl;
    }

    for (auto t : coords){
        Eigen::Vector3d curr(t[0], t[1], t[2]);
        coordinates.push_back(curr);
    }

    for (auto t : orients){
        Rotation rot(t[0], t[1], t[2]);
        orientations.push_back(rot);
    }

    Eigen::Vector3d c(cell_[0], cell_[1], cell_[2]);
    cell = c;

    N = coordinates.size();

    // automatically set all types to 0
    types.resize(N);
    std::fill(types.begin(), types.end(), 0);

    // by default set time_stamp to 0
    time_stamp = 0;
}

// getters
int Frame::get_N() { return N; }
int Frame::get_time_stamp() { return time_stamp; }
std::vector<Eigen::Vector3d> Frame::get_coordinates() { return coordinates; }
std::vector<std::vector<double> > Frame::get_coordinates_as_list(){
    std::vector<std::vector<double> > res;
    for(int i = 0; i < N; i++){
        std::vector<double> row;
        for(int j = 0; j < 3; j++){
            row.push_back(coordinates[i](j));
        }
        res.push_back(row);
    }
    return res;
}
std::vector<Rotation> Frame::get_orientations() { return orientations; }
std::vector<std::vector<double> > Frame::get_orientations_as_list_of_angles(){
    std::vector<std::vector<double> > res;
    for(int i = 0; i < N; i++){
        res.push_back(orientations[i].get_zxz_ccw_angles());
    }
    return res;
}
Eigen::Vector3d Frame::get_cell() { return cell; }
std::vector<double> Frame::get_cell_as_list(){
    std::vector<double> res;
    for(int i = 0; i < 3; i++) res.push_back(cell(i));
    return res;
}
std::vector<int> Frame::get_types() { return types; }
std::vector<std::vector<int> > Frame::get_bond_list(){ return bond_list; }

// setters
void Frame::set_time_stamp(int fn) { time_stamp = fn; }
void Frame::set_types(std::vector<int> ty) { types = ty; }

// write xyz file
void Frame::write_xyz(std::string file_name){
    std::ofstream out_file;
    out_file.open(file_name);
    out_file << N << std::endl;
    out_file << "time = " << time_stamp;
    out_file << ", Lx= " << cell(0);
    out_file << ", Ly= " << cell(1);
    out_file << ", L!= " << cell(2) << std::endl;

    for (int i = 0; i < N; i++){
        // TODO: switch to unordered_map
        if(types[i] == 0) out_file << "N ";
        else if (types[i] == 1) out_file << "C ";
        else if (types[i] == 2) out_file << "O ";
        else if (types[i] == 3) out_file << "H ";
        else out_file << "S ";

        out_file << coordinates[i](0) << " " << coordinates[i](1) << " " << coordinates[i](2) << std::endl;
    }
    out_file.close();
}

// interactions
Eigen::Vector3d Frame::get_ij_displacement(int i, int j){
    Eigen::Vector3d disp = coordinates[i] - coordinates[j];

    // apply pbc
    Eigen::Vector3d h_cell = cell / 2;
    for (int i = 0; i < 3; i++){
        while (disp(i) > h_cell(i)) disp(i) -= cell(i);
        while (disp(i) < -h_cell(i)) disp(i) += cell(i);
    }

    return disp;
}
bool Frame::does_ij_interact(int i, int j, Patches& p){
    // check dist
    double max_lambda = p.get_max_lambda();
    Eigen::Vector3d disp = this->get_ij_displacement(i, j);
    if (disp.norm() > max_lambda) return false;
    if (disp.norm() < 1) return false; // overlap

    // check orientations
    disp.normalize();
    double omega_i, omega_j;
    for(int pi = 0; pi < p.get_n_patch(); pi++){
        if(p.get_type(pi) != types[i]) continue;
        Eigen::Vector3d pi_v = p.get_vector(pi);
        pi_v = orientations[i].rotate_vec(pi_v);
        for(int pj = 0; pj < p.get_n_patch(); pj++){
            if(p.get_type(pj) != types[j]) continue;
            Eigen::Vector3d pj_v = p.get_vector(pj);
            pj_v = orientations[i].rotate_vec(pj_v);

            omega_i = pi_v.transpose() * disp;
            omega_j = pj_v.transpose() * (-1*disp);

            if (omega_i >= p.get_cos_delta(pi) && omega_j >= p.get_cos_delta(pj)) return true;
        }
    }
    return false;
}
void Frame::determine_bond_list(Patches& p){
    bond_list.clear();
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            if(this->does_ij_interact(i, j, p)){
                std::vector<int> pair = {i, j};
                bond_list.push_back(pair);
            }
        }
    }
}