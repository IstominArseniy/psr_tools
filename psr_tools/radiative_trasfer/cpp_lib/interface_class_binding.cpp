#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "CppInterfaceClass.h"

namespace py = pybind11;

PYBIND11_MODULE(loki_python_binding, m, py::mod_gil_not_used()) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    py::class_<CppInterface>(m, "ProfileCalculator")
        .def(py::init<const std::map<std::string, double>, const std::map<std::string, double> >())
        .def("find_ILVPA", &CppInterface::find_ILVPA)
        .def("find_I", &CppInterface::find_I)
        .def("get_rho", &CppInterface::get_rho)
        .def("calculate_profile", &CppInterface::calculate_profile);

}