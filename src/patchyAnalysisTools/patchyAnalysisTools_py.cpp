#include "frame.hpp"
#include "trajectory.hpp"
#include "patches.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

void init_patches(py::module_ &m){
    py::class_<Patches>(m, "Patches", py::dynamic_attr())
        .def(py::init<>())
        .def("add_patch", py::overload_cast<double, double, double, std::vector<double>, int>(&Patches::add_patch))
        .def("read_from_file", &Patches::read_from_file)
        .def("make_ij_interact", &Patches::make_ij_interact)
        .def("make_all_patches_adjacent", &Patches::make_all_patches_adjacent)
        .def("are_ij_adjacent", &Patches::are_ij_adjacent)
        ;
}

PYBIND11_MODULE(_patchyAnalysisTools, m){
    init_frame(m);
    init_trajectory(m);
}