/*
C++ module with pybind11 to analyze patchy particle simulation data.
*/

#include <pybind11/pybind11.h>
#include "include/trajectory.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_trajectory(py::module_ &);
void init_frame(py::module_ &);

PYBIND11_MODULE(patchyAnalysisTools_cpp, m){
    init_trajectory(m);
    init_frame(m);

    #ifdef VERSION_INFO
        m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
        m.attr("__version__") = "dev";
    #endif
}