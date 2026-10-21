#include "beam.hpp"
#include <iostream>
#include <limits>

int main() {
    const auto e = bucklemap::element(2.0, 3.0);
    const std::array<double, 4> translation{1, 0, 1, 0};
    const std::array<double, 4> rotation{0, 1, 2, 1};
    for (int i = 0; i < 4; ++i) {
        double kt = 0, kr = 0, gt = 0;
        for (int j = 0; j < 4; ++j) {
            if (e.elastic[4 * i + j] != e.elastic[4 * j + i] ||
                e.geometric[4 * i + j] != e.geometric[4 * j + i])
                return 1;
            kt += e.elastic[4 * i + j] * translation[j];
            kr += e.elastic[4 * i + j] * rotation[j];
            gt += e.geometric[4 * i + j] * translation[j];
        }
        if (std::abs(kt) + std::abs(kr) + std::abs(gt) > 1e-13)
            return 2;
    }
    for (double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        try {
            bucklemap::element(bad, 1.0);
            return 3;
        } catch (const std::invalid_argument &) {
        }
    }
    std::cout << "Symmetry, rigid modes and invalid input checks passed\n";
    return 0;
}
