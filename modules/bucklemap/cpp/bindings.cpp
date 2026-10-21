#include "beam.hpp"
#include <algorithm>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
PYBIND11_MODULE(_core, m) {
    m.doc() = "Exactly integrated cubic Hermite beam matrices (C++17)";
    m.def(
        "element_matrices",
        [](double length, double rigidity) {
            const auto e = bucklemap::element(length, rigidity);
            py::array_t<double> k({4, 4}), g({4, 4});
            std::copy(e.elastic.begin(), e.elastic.end(), k.mutable_data());
            std::copy(e.geometric.begin(), e.geometric.end(), g.mutable_data());
            return py::make_tuple(k, g);
        },
        py::arg("length"), py::arg("rigidity"));
}
