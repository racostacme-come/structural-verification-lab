#include "assembly.h"
#include <algorithm>
#include <cmath>
#include <limits>

int tm_abi_version(void) { return 1; }

int tm_assemble(size_t n, const double* x, const double* young,
                const double* area, const double* alpha, const double* t,
                double* diagonal, double* off_diagonal, double* load) {
    if (n == 0 || n == std::numeric_limits<size_t>::max() || !x || !young ||
        !area || !alpha || !t || !diagonal || !off_diagonal || !load) return 1;
    for (size_t i = 0; i <= n; ++i)
        if (!std::isfinite(x[i])) return 2;
    for (size_t e = 0; e < n; ++e) {
        if (!(x[e+1] > x[e]) || !std::isfinite(young[e]) || young[e] <= 0 ||
            !std::isfinite(area[e]) || area[e] <= 0 ||
            !std::isfinite(alpha[e]) || !std::isfinite(t[e])) return 2;
    }
    std::fill(diagonal, diagonal+n+1, 0.0);
    std::fill(load, load+n+1, 0.0);
    for (size_t e = 0; e < n; ++e) {
        const double h = x[e+1]-x[e];
        const double ea = young[e]*area[e];
        const double k = ea/h;
        const double f = ea*alpha[e]*t[e];
        if (!std::isfinite(h) || !std::isfinite(ea) || !std::isfinite(k) ||
            !std::isfinite(f) || k <= 0) return 3;
        diagonal[e] += k;
        diagonal[e+1] += k;
        off_diagonal[e] = -k;
        load[e] -= f;
        load[e+1] += f;
        if (!std::isfinite(diagonal[e]) || !std::isfinite(diagonal[e+1]) ||
            !std::isfinite(load[e]) || !std::isfinite(load[e+1])) return 3;
    }
    return 0;
}
