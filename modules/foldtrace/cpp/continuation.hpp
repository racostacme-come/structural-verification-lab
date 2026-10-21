#pragma once

#include <cmath>
#include <stdexcept>

namespace foldtrace {

// Symmetric two-bar truss, linear engineering-strain axial law, exact geometry.
// r=h/a, q=w/h, load=P/(2 EA h/L0), energy=U/(EA a).
class Truss {
  public:
    explicit Truss(double ratio) : r_(ratio), l0_(std::hypot(1.0, ratio)) {
        if (!std::isfinite(ratio) || ratio < 0.01 || ratio > 10.0)
            throw std::invalid_argument("height/span ratio must be in [0.01, 10]");
    }
    double ratio() const { return r_; }
    double force(double q) const {
        const double l = length(q);
        // Rationalization avoids subtracting near-equal lengths at q=0 and 2.
        return (1.0 - q) * r_ * r_ * q * (2.0 - q) / (l * (l0_ + l));
    }
    double tangent(double q) const {
        const double l = length(q);
        return 1.0 - l0_ / (l * l * l);
    }
    double energy(double q) const {
        const double l = length(q);
        const double dl = r_ * r_ * q * (q - 2.0) / (l + l0_);
        return dl * dl / l0_;
    }

  private:
    double r_, l0_;
    double length(double q) const {
        if (!std::isfinite(q) || q < -10.0 || q > 10.0)
            throw std::invalid_argument("displacement q must be finite and in [-10, 10]");
        return std::hypot(1.0, r_ * (1.0 - q));
    }
};

struct Correction {
    double q, load, residual, constraint;
    int iterations;
    bool converged;
};

// Orthogonal hyperplane corrector in z=(q, load/load_scale).
// Predictor and tangent are held fixed during Newton iterations.
inline Correction correct(const Truss &model, double qp, double lp, double tq, double tl,
                          double scale, double tolerance, int max_iterations) {
    if (!std::isfinite(qp) || !std::isfinite(lp) || !std::isfinite(tq) || !std::isfinite(tl) ||
        !std::isfinite(scale) || scale <= 0.0 || !std::isfinite(tolerance) ||
        tolerance <= 0.0 || max_iterations < 1 || std::abs(std::hypot(tq, tl) - 1.0) > 1e-10)
        throw std::invalid_argument("invalid corrector settings or non-unit tangent");
    double q = qp, load = lp;
    double residual = 0.0, constraint = 0.0;
    for (int iteration = 0; iteration <= max_iterations; ++iteration) {
        if (!std::isfinite(q) || std::abs(q) > 10.0 || !std::isfinite(load))
            return {q, load, INFINITY, INFINITY, iteration, false};
        residual = (model.force(q) - load) / scale;
        constraint = tq * (q - qp) + tl * (load - lp) / scale;
        if (!std::isfinite(residual) || !std::isfinite(constraint))
            return {q, load, residual, constraint, iteration, false};
        if (std::hypot(residual, constraint) <= tolerance)
            return {q, load, residual, constraint, iteration, true};
        if (iteration == max_iterations)
            break;
        const double k = model.tangent(q) / scale;
        const double det = k * tl + tq;
        if (!std::isfinite(det) || std::abs(det) < 1e-14)
            return {q, load, residual, constraint, iteration, false};
        const double dq = (-tl * residual - constraint) / det;
        const double dz = (tq * residual - k * constraint) / det;
        q += dq;
        load += scale * dz;
    }
    return {q, load, residual, constraint, max_iterations, false};
}

} // namespace foldtrace
