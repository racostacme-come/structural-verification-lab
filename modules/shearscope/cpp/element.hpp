#pragma once
#include <array>
#include <cmath>
#include <stdexcept>

namespace shearscope {
using Matrix = std::array<double, 16>;
struct Element {
    Matrix stiffness{};
    Matrix mass{};
};

// DOF order: transverse displacement and rotation at each of two nodes.
inline Element element(double length, double ei, double shear, double rho_a, double rho_i,
                       bool reduced) {
    for (double value : {length, ei, shear, rho_a, rho_i}) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw std::invalid_argument("element parameters must be finite and positive");
        }
    }
    Element out;
    const std::array<double, 4> curvature{0.0, -1.0 / length, 0.0, 1.0 / length};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            out.stiffness[4 * i + j] = ei * length * curvature[i] * curvature[j];
        }
    }
    const double point = 1.0 / std::sqrt(3.0);
    const int count = reduced ? 1 : 2;
    for (int q = 0; q < count; ++q) {
        const double xi = reduced ? 0.0 : (q == 0 ? -point : point);
        const double weight = reduced ? length : length / 2.0;
        const std::array<double, 4> gamma{-1.0 / length, -(1.0 - xi) / 2.0, 1.0 / length,
                                          -(1.0 + xi) / 2.0};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                out.stiffness[4 * i + j] += shear * weight * gamma[i] * gamma[j];
            }
        }
    }
    // Consistent mass is always integrated exactly, including rotary inertia.
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            const double factor = length * (a == b ? 2.0 : 1.0) / 6.0;
            out.mass[4 * (2 * a) + 2 * b] = rho_a * factor;
            out.mass[4 * (2 * a + 1) + 2 * b + 1] = rho_i * factor;
        }
    }
    for (const auto &matrix : {out.stiffness, out.mass}) {
        for (double value : matrix) {
            if (!std::isfinite(value)) {
                throw std::overflow_error("element matrix overflow; rescale physical units");
            }
        }
    }
    return out;
}
} // namespace shearscope
