#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "include/trajectory.hpp"

namespace py = pybind11;

Trajectory::Trajectory(std::vector<Frame> list_of_frames_){ list_of_frames = list_of_frames_; }

Frame Trajectory::get_frame(int frame_idx){ return list_of_frames[frame_idx];}

void init_trajectory(py::module_ &m){
    py::class_<Trajectory>(m, "Trajectory", py::dynamic_attr())
        .def(py::init<std::vector<Frame>>())
        .def("get_frame", &Trajectory::get_frame)
        ;
}