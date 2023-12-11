/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#pragma once
#include "frame.hpp"

class Trajectory{
    std::vector<Frame> frames;
    Patches patches;
public:
    // constructor
    Trajectory();
    void add_frame(const Frame& frame);

    // setters
    void set_patches(Patches& patches_);

    // getters
    Frame get_frame(int idx);
};