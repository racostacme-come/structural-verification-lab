#include "element.hpp"
#include <iostream>

int main() {
    for (bool reduced : {false, true}) {
        const auto result = shearscope::element(2.0, 3.0, 7.0, 5.0, 0.2, reduced);
        for (const auto &mode :
             {std::array<double, 4>{1, 0, 1, 0}, std::array<double, 4>{0, 1, 2, 1}}) {
            for (int i = 0; i < 4; ++i) {
                double force = 0;
                for (int j = 0; j < 4; ++j) {
                    force += result.stiffness[4 * i + j] * mode[j];
                    if (std::abs(result.stiffness[4 * i + j] - result.stiffness[4 * j + i]) >
                        1e-12) {
                        return 1;
                    }
                }
                if (std::abs(force) > 1e-12) {
                    return 2;
                }
            }
        }
        if (std::abs(result.mass[0] + result.mass[2] + result.mass[8] + result.mass[10] - 10.0) >
            1e-12) {
            return 3;
        }
    }
    try {
        shearscope::element(-1, 1, 1, 1, 1, true);
        return 4;
    } catch (const std::invalid_argument &) {
    }
    std::cout << "Rigid motions, symmetry, total mass and invalid input verified\n";
    return 0;
}
