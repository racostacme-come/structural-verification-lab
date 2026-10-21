#pragma once
#include <array>
#include <cmath>
#include <stdexcept>

namespace bucklemap {
using Matrix = std::array<double, 16>;
struct Element {
    Matrix elastic;
    Matrix geometric;
};

// DOF order: transverse displacement and slope at each end.
// geometric is positive for a UNIT compressive dead load: (K - P G) q = 0.
inline Element element(double length, double rigidity) {
    if (!std::isfinite(length) || length <= 0 || !std::isfinite(rigidity) || rigidity <= 0)
        throw std::invalid_argument("length and EI must be finite and positive");
    const double h = length;
    const double h2 = h * h;
    Element e{{12, 6 * h, -12, 6 * h, 6 * h, 4 * h2, -6 * h, 2 * h2, -12, -6 * h, 12, -6 * h, 6 * h,
               2 * h2, -6 * h, 4 * h2},
              {36, 3 * h, -36, 3 * h, 3 * h, 4 * h2, -3 * h, -h2, -36, -3 * h, 36, -3 * h, 3 * h,
               -h2, -3 * h, 4 * h2}};
    for (auto &v : e.elastic) {
        v *= rigidity / (h * h2);
        if (!std::isfinite(v))
            throw std::overflow_error("elastic matrix exceeds floating-point range");
    }
    for (auto &v : e.geometric) {
        v /= 30 * h;
        if (!std::isfinite(v))
            throw std::overflow_error("geometric matrix exceeds floating-point range");
    }
    if (e.elastic[0] <= 0 || e.geometric[0] <= 0)
        throw std::overflow_error("element scaling underflows floating-point range");
    return e;
}
} // namespace bucklemap
