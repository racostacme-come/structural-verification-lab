#pragma once
#include <cmath>
#include <stdexcept>

namespace yieldledger {
struct State {
    double plastic_strain = 0.0;
    double accumulated = 0.0;
};

struct Response {
    double stress;
    double tangent;
    double plastic_increment;
    State state;
};

class Material {
  public:
    const double elastic, yield_stress, isotropic, kinematic;

    Material(double e, double y, double hi, double hk)
        : elastic(e), yield_stress(y), isotropic(hi), kinematic(hk) {
        if (!std::isfinite(e) || !std::isfinite(y) || !std::isfinite(hi) ||
            !std::isfinite(hk) || e <= 0 || y <= 0 || hi < 0 || hk < 0 ||
            !std::isfinite(e + hi + hk)) {
            throw std::invalid_argument("Require finite E,y > 0 and Hi,Hk >= 0.");
        }
    }

    static void validate(double strain, const State &s) {
        if (!std::isfinite(strain) || !std::isfinite(s.plastic_strain) ||
            !std::isfinite(s.accumulated) || s.accumulated < 0 ||
            std::abs(s.plastic_strain) > s.accumulated + 1e-14) {
            throw std::invalid_argument("Require finite strain and alpha >= |plastic strain|.");
        }
    }

    Response update(double strain, const State &old = {}) const {
        validate(strain, old);
        const double trial = elastic * (strain - old.plastic_strain);
        const double relative = trial - kinematic * old.plastic_strain;
        const double f = std::abs(relative) - yield_stress - isotropic * old.accumulated;
        if (!std::isfinite(trial) || !std::isfinite(relative) || !std::isfinite(f)) {
            throw std::overflow_error("Constitutive trial overflow.");
        }
        if (f <= 0) {
            return {trial, elastic, 0.0, old};
        }
        const double dg = f / (elastic + isotropic + kinematic);
        const double sign = std::copysign(1.0, relative);
        State next{old.plastic_strain + sign * dg, old.accumulated + dg};
        const double stress = trial - elastic * sign * dg;
        const double tangent =
            elastic * ((isotropic + kinematic) / (elastic + isotropic + kinematic));
        if (!std::isfinite(stress) || !std::isfinite(next.accumulated) ||
            !std::isfinite(next.plastic_strain)) {
            throw std::overflow_error("Constitutive update overflow.");
        }
        return {stress, tangent, dg, next};
    }

    double energy(double strain, const State &s) const {
        validate(strain, s);
        const double ee = strain - s.plastic_strain;
        const double value = 0.5 * elastic * ee * ee +
                             0.5 * isotropic * s.accumulated * s.accumulated +
                             0.5 * kinematic * s.plastic_strain * s.plastic_strain;
        if (!std::isfinite(value)) {
            throw std::overflow_error("Stored energy overflow.");
        }
        return value;
    }
};
} // namespace yieldledger
