#include "continuation.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;
PYBIND11_MODULE(_core, m) {
    m.doc() = "Original C++17 truss mechanics and augmented Newton corrector";
    py::class_<foldtrace::Truss>(m, "Truss")
        .def(py::init<double>(), py::arg("ratio") = 0.2)
        .def_property_readonly("ratio", &foldtrace::Truss::ratio)
        .def("force", &foldtrace::Truss::force, py::arg("q"))
        .def("tangent", &foldtrace::Truss::tangent, py::arg("q"))
        .def("energy", &foldtrace::Truss::energy, py::arg("q"));
    py::class_<foldtrace::Correction>(m, "Correction")
        .def_readonly("q", &foldtrace::Correction::q)
        .def_readonly("load", &foldtrace::Correction::load)
        .def_readonly("residual", &foldtrace::Correction::residual)
        .def_readonly("constraint", &foldtrace::Correction::constraint)
        .def_readonly("iterations", &foldtrace::Correction::iterations)
        .def_readonly("converged", &foldtrace::Correction::converged);
    m.def("correct", &foldtrace::correct, py::arg("model"), py::arg("qp"), py::arg("lp"),
          py::arg("tq"), py::arg("tl"), py::arg("scale"), py::arg("tolerance"),
          py::arg("max_iterations"), py::call_guard<py::gil_scoped_release>());
}
