#include "trajectory.hpp"

Trajectory::Trajectory(std::vector<Frame> list_of_frames_){ list_of_frames = list_of_frames_; }

Frame Trajectory::get_frame(int frame_idx){ return list_of_frames[frame_idx];}