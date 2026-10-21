#include "element.hpp"
#include <algorithm>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
PYBIND11_MODULE(_core, module) {
    module.doc() = "Original two-node Timoshenko element kernel";
    module.def(
        "element",
        [](double length, double ei, double shear, double rho_a, double rho_i, bool reduced) {
            const auto result = shearscope::element(length, ei, shear, rho_a, rho_i, reduced);
            py::array_t<double> stiffness({4, 4});
            py::array_t<double> mass({4, 4});
            std::copy(result.stiffness.begin(), result.stiffness.end(), stiffness.mutable_data());
            std::copy(result.mass.begin(), result.mass.end(), mass.mutable_data());
            return py::make_tuple(stiffness, mass);
        },
        py::arg("length"), py::arg("ei"), py::arg("shear"), py::arg("rho_a"), py::arg("rho_i"),
        py::arg("reduced") = true);
}
