#include "material.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace yieldledger;

PYBIND11_MODULE(_core, m) {
    m.doc() = "Immutable uniaxial mixed-hardening return mapping (C++17).";
    py::class_<State>(m, "State")
        .def(py::init<double, double>(), py::arg("plastic_strain") = 0.0,
             py::arg("accumulated") = 0.0)
        .def_readonly("plastic_strain", &State::plastic_strain)
        .def_readonly("accumulated", &State::accumulated);
    py::class_<Response>(m, "Response")
        .def_readonly("stress", &Response::stress)
        .def_readonly("tangent", &Response::tangent)
        .def_readonly("plastic_increment", &Response::plastic_increment)
        .def_readonly("state", &Response::state);
    py::class_<Material>(m, "Material")
        .def(py::init<double, double, double, double>(), py::arg("elastic"),
             py::arg("yield_stress"), py::arg("isotropic") = 0.0, py::arg("kinematic") = 0.0)
        .def_readonly("elastic", &Material::elastic)
        .def_readonly("yield_stress", &Material::yield_stress)
        .def_readonly("isotropic", &Material::isotropic)
        .def_readonly("kinematic", &Material::kinematic)
        .def("update", &Material::update, py::arg("strain"), py::arg("old") = State{})
        .def("energy", &Material::energy, py::arg("strain"), py::arg("state"));
}
