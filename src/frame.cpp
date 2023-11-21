
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "include/frame.hpp"

namespace py = pybind11;

// constructors
Frame::Frame(std::vector<Vec3> coordinates_, std::vector<Rotation> orientations_){
    if (coordinates.size() != orientations.size()) throw "coordinates and orientations have different sizes.";
    coordinates = coordinates_;
    orientations = orientations_;

    N = coordinates.size();
    std::fill(types.begin(), types.end(), 0);
    frame_num = 0;
}

Frame::Frame(std::vector<std::vector<double> > coords, std::vector<std::vector<double> > orients){
    if (coords.size() != orients.size()) throw "coordinates and orientations have different sizes.";
    N = coords.size();
    for( auto t : coords){
        Vec3 curr(t);
        coordinates.push_back(curr);
    }
    for(auto t : orients){
        Rotation rot(t[1], t[2], t[3]);
        orientations.push_back(rot);
    }

    std::fill(types.begin(), types.end(), 0);
    frame_num = 0;
}

// accessing private props (do we even want these?)
void Frame::set_coordinates(std::vector<Vec3> coords){ 
    if( orientations.size() > 0 && orientations.size() != coords.size()) throw "coordinates size is different from orientations"; 

    coordinates = coords;
}
void Frame::set_coordinates(std::vector<std::vector<double> > coords){ 
    if( orientations.size() > 0 && orientations.size() != coords.size()) throw "coordinates size is different from orientations"; 
    coordinates.clear();
    for (auto t : coords){
        Vec3 curr(t);
        coordinates.push_back(curr);
    }
}

int Frame::get_N(){ return N; }
void Frame::set_frame_num(int num){frame_num = num;}
int Frame::get_frame_num(){return frame_num;}

std::vector<Vec3> Frame::get_coordinates(){ return coordinates; }
std::vector<std::vector<double> > Frame::get_orientations_as_vector_of_angles(){
    std::vector<std::vector<double> > res;
    for( auto& t : orientations){
        res.push_back(t.get_zxz_ccw_angles());
    }
    return res;
}

// again, do we event want these?
void Frame::set_orientations(std::vector<Rotation> orient){
    if(coordinates.size() > 0 && coordinates.size() != orient.size()) throw "orientations size is different from coordinates";
    orientations = orient;
}
void Frame::set_orientations(std::vector<Vec3> orient){
    if(coordinates.size() > 0 && coordinates.size() != orient.size()) throw "orientations size is different from coordinates";
    std::vector<double> avec;
    for(auto t : orient){
        avec = t.get_vec();
        Rotation curr(avec[0], avec[2], avec[3]);
        orientations.push_back(curr);
    }
}
std::vector<Rotation> Frame::get_orientations() { return orientations; }

void Frame::write_xyz(std::string file_name){
    std::ofstream out_file;
    std::vector<double> c = cell.get_vec();
    out_file.open(file_name);
    out_file << N << std::endl;
    out_file << "frame = " << frame_num;
    out_file << ", Lx= " << c[0];
    out_file << ", Ly= " << c[1];
    out_file << ", Lz= " << c[2] << std::endl;

    for(int i = 0; i < N; i ++){
        // TODO: switch to unordered_map
        if(types[i] == 0) out_file << "N ";
        else if (types[i] == 1) out_file << "C ";
        else if (types[i] == 2) out_file << "O ";
        else if (types[i] == 3) out_file << "H ";
        else out_file << "S ";

        std::vector<double> vec = coordinates[i].get_vec();
        out_file << vec[0] << " " << vec[1] << " " << vec[2] << " " << std::endl;
    }
    out_file.close();
}

void Frame::set_patch_info(Patches patches_) { patches = patches;}
void Frame::set_patch_info(std::string file_name){
    patches.read_from_file(file_name);
}

bool Frame::does_i_j_interact(int i, int j, double tol){
    // check dist
    double max_lambda = this->patches.get_max_lambda();
    Vec3 disp = this->get_i_j_displacement(i, j);
    if( disp.get_norm_sq() > max_lambda * max_lambda) return false;

    // check orientations
    disp.normalize();
    Vec3 disp_flipped = disp.get_flipped();
    double omega_i, omega_j;
    for(int pi; pi < patches.get_n_patch(); pi++){
        if(!(patches.check_type(types[i], pi))) continue;
        Vec3 pi_v = patches.get_vector(pi);
        pi_v.rotate(orientations[i]);
        for(int pj; pj < patches.get_n_patch(); pj++){
            if(!(patches.check_type(types[j], pj))) continue;
            if(!patches.are_ij_adjacent(pi, pj)) continue;
            Vec3 pj_v = patches.get_vector(pj);
            pj_v.rotate(orientations[j]);
            omega_i = pi_v.dot(disp);
            omega_j = pj_v.dot(disp_flipped);

            if(omega_i+tol >= patches.get_cos_delta(pi) && omega_j+tol >= patches.get_cos_delta(pj)) return true;
        }
    }
    return false;
}

Vec3 Frame::get_i_j_displacement(int i, int j){
    Vec3 disp = coordinates[i].diff(coordinates[j]);
    // apply periodic boundary conditions
    disp.apply_pbc(cell);
    return disp;
}

std::vector<std::vector<int> > Frame::get_list_of_interacting_pairs(double tol){
    if(patches.get_n_patch() == 0) throw "no patches to calculate energy with.";

    bonds_calculated.clear();
    std::fill(bonds_calculated.begin(), bonds_calculated.end(), 0);

    std::vector<std::vector<int> > res;

    // no neighbor lists because this is not meant to be repeated a ton
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            if(does_i_j_interact(i, j)){
                bonds_calculated[i] += 1;
                bonds_calculated[j] += 1;
                std::vector<int> pair;
                pair.push_back(i);
                pair.push_back(j);
                res.push_back(pair);
            }
        }
    }
    return res;
}

bool Frame::check_calculated_bonds_against_bond_numbers(){
    if(bonds_calculated.size() == 0) throw "bonds have not been calculated";
    for(int i = 0; i < N; i++){
        if (bonding[i] != bonds_calculated[i]) return false;
    }
    return true;
}

void Frame::set_cluster_info(){
    std::vector<std::vector<int> > interacting_pairs = this->get_list_of_interacting_pairs();
    cluster_info.set_bonds(interacting_pairs);
    cluster_info.clusters_from_interacting_pairs(N);

}

void init_frame(py::module_ &m){
    py::class_<Frame>(m, "Frame", py::dynamic_attr())
        .def(py::init<std::vector<std::vector<double> >, std::vector<std::vector<double> > >())
        ;
}