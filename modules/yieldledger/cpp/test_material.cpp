#include "material.hpp"
#include <iostream>
#include <limits>

void check(bool condition) {
    if (!condition) {
        throw std::runtime_error("material invariant failed");
    }
}

int main() {
    using namespace yieldledger;
    Material m(200000, 250, 2000, 3000);
    const auto elastic = m.update(0.0001);
    check(std::abs(elastic.stress - 20) < 1e-12);
    auto s = State{};
    double previous_strain = 0;
    for (int i = 0; i < 2000; ++i) {
        const double strain = 0.008 * std::sin(0.03 * i);
        const auto r = m.update(strain, s);
        const double dep = r.state.plastic_strain - s.plastic_strain;
        const double da = r.state.accumulated - s.accumulated;
        const double dee = strain - previous_strain - dep;
        const double defect = 0.5 * m.elastic * dee * dee + 0.5 * m.isotropic * da * da +
                              0.5 * m.kinematic * dep * dep;
        const double balance = r.stress * (strain - previous_strain) -
                               (m.energy(strain, r.state) - m.energy(previous_strain, s)) -
                               m.yield_stress * da - defect;
        check(std::abs(balance) < 1e-10);
        check(da >= 0);
        check(std::abs(r.stress - m.kinematic * r.state.plastic_strain) - m.yield_stress -
                  m.isotropic * r.state.accumulated <
              1e-9);
        s = r.state;
        previous_strain = strain;
    }
    bool rejected = false;
    try {
        m.update(std::numeric_limits<double>::quiet_NaN());
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected);
    std::cout
        << "2000 cyclic updates: energy identity, yield admissibility, irreversibility OK\n";
}
