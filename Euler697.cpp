#include <boost/math/special_functions/gamma.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {

long double log10_c_for_n(const long double n) {
    const long double ln_c = boost::math::gamma_p_inv(n, 0.75L);
    return ln_c / std::log(10.0L);
}

}  // namespace

int main() {
    const long double ln_c_100 = boost::math::gamma_p_inv(100.0L, 0.75L);
    const long double log10_c_100 = ln_c_100 / std::log(10.0L);
    assert(std::fabsl(log10_c_100 - 46.27L) < 0.01L);
    assert(std::fabsl(boost::math::gamma_q(100.0L, ln_c_100) - 0.25L) < 1e-15L);

    const long double n = 10'000'000.0L;
    const long double ln_c = boost::math::gamma_p_inv(n, 0.75L);
    assert(std::fabsl(boost::math::gamma_q(n, ln_c) - 0.25L) < 1e-12L);

    const long double answer = ln_c / std::log(10.0L);
    std::cout << std::fixed << std::setprecision(2) << static_cast<double>(answer) << '\n';
    return 0;
}

