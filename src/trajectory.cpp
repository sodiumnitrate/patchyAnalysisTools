/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#include "trajectory.hpp"

// constructor
Trajectory::Trajectory(){}

// add frame
void Trajectory::add_frame(const Frame& frame){
    frames.push_back(frame);
}

// setters
void Trajectory::set_patches(Patches& patches_){
    patches = patches_;
}

// getters
Frame Trajectory::get_frame(int idx){
    return frames[idx];
}