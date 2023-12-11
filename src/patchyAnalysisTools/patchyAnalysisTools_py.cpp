/*
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.
*/

#include "frame.hpp"
#include "patches.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_frame(py::module_ &m){
    py::class_<Frame>(m, "Frame")
        .def(py::init<std::vector<std::vector<double> >, std::vector<std::vector<double> >, std::vector<double> >())
        .def("get_N", &Frame::get_N)
        .def("get_time_stamp", &Frame::get_time_stamp)
        .def("get_coordinates", &Frame::get_coordinates_as_list)
        .def("get_orientations_as_angles", &Frame::get_orientations_as_list_of_angles)
        .def("get_types", &Frame::get_types)
        .def("get_bond_list", &Frame::get_bond_list)
        .def("set_time_stamp", &Frame::set_time_stamp)
        .def("set_types", &Frame::set_types)
        .def("write_xyz", &Frame::write_xyz)
        .def("does_ij_interact", &Frame::does_ij_interact)
        .def("determine_bond_list", &Frame::determine_bond_list)
        .def("determine_clusters", &Frame::determine_clusters)
        .def("determine_percolation", &Frame::determine_percolation)
        .def("get_clusters", &Frame::get_clusters)
        ;
}

void init_patches(py::module_ &m){
    py::class_<Patches>(m, "Patches")
        .def(py::init<>())
        .def("add_patch", py::overload_cast<double, double, double, std::vector<double>, int>(&Patches::add_patch))
        .def("read_from_file", &Patches::read_from_file)
        .def("get_n_patch", &Patches::get_n_patch)
        .def("get_max_lambda", &Patches::get_max_lambda)
        .def("make_all_patches_adjacent", &Patches::make_all_patches_adjacent)
        ;
}

void init_clusters(py::module_ &m){
    py::class_<Clusters>(m, "Clusters")
        .def(py::init<>())
        .def("get_clusters", &Clusters::get_clusters)
        .def("get_bond_list", &Clusters::get_bond_list)
        ;
}

PYBIND11_MODULE(_patchyAnalysisTools, m){
    init_frame(m);
    init_patches(m);
    init_clusters(m);
}