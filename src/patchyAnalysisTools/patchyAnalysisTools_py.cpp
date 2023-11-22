#include "frame.hpp"
#include "trajectory.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_frame(py::module_ &m){
    py::class_<Frame>(m, "Frame", py::dynamic_attr())
        .def(py::init<std::vector<std::vector<double> >, std::vector<std::vector<double> > , std::vector<double> >())
        .def("get_N", &Frame::get_N)
        .def("get_frame_num", &Frame::get_frame_num)
        .def("set_frame_num", &Frame::set_frame_num)
        .def_property("frame_num", &Frame::get_frame_num, &Frame::set_frame_num)
        .def("get_coordinates", &Frame::get_coordinates_as_list)
        .def("get_orientations", &Frame::get_orientations_as_vector_of_angles)
        .def("write_xyz", &Frame::write_xyz)
        .def("get_cell", &Frame::get_cell_as_list)
        .def("get_displacement", &Frame::get_i_j_displacement_as_vec)
        ;
}

void init_trajectory(py::module_ &m){
    py::class_<Trajectory>(m, "Trajectory", py::dynamic_attr())
        .def(py::init<std::vector<Frame>>())
        .def("get_frame", &Trajectory::get_frame)
        ;
}

PYBIND11_MODULE(_patchyAnalysisTools, m){
    init_frame(m);
    init_trajectory(m);
}