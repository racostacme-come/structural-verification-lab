#include "continuation.hpp"
#include <iostream>
#include <limits>

int main() {
    int failures = 0;
    const auto check = [&failures](bool ok, const char *message) {
        if (!ok) {
            std::cerr << message << '\n';
            ++failures;
        }
    };
    const foldtrace::Truss model(0.2);
    check(model.force(0.0) == 0.0 && model.force(2.0) == 0.0, "stress-free states");
    check(model.tangent(1.0) < 0.0, "negative stiffness at flattened state");
    const double eps = 1e-5;
    check(std::abs((model.force(0.3 + eps) - model.force(0.3 - eps)) / (2 * eps) -
                   model.tangent(0.3)) < 1e-10,
          "consistent tangent");
    const double k = model.tangent(0.0) / 0.04;
    const double tq = 1.0 / std::hypot(1.0, k), tl = k * tq;
    const auto result =
        foldtrace::correct(model, 0.1 * tq, 0.004 * tl, tq, tl, 0.04, 1e-12, 12);
    check(result.converged && std::abs(model.force(result.q) - result.load) < 1e-12,
          "augmented Newton equilibrium");
    bool rejected = false;
    try {
        model.force(std::numeric_limits<double>::quiet_NaN());
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "reject non-finite input");
    return failures == 0 ? 0 : 1;
}
