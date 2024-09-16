#pragma once
#include "frame.hpp"

class Trajectory{
    std::vector<Frame> list_of_frames;
    int N_frames;
    // TODO: where to hold the patch info
public:
    Trajectory(std::vector<Frame> list_of_frames_);
    //Trajectory(std::string file_name);
    Frame get_frame(int frame_idx);
    //void write_xyz(std::string file_name);
    //void set_patch_info(Patches patches_);
    //void set_patch_info(std::string file_name);
    //void set_energy();
    //void set_energy(std::vector<double> energy);
    //std::vector<double> get_energy();

    /*
    TODO:
    - get_en_autocorr
    - position_autocorr
    - read_trajectory
    - get_density_distribution
    */
};