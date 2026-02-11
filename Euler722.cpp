#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace {

using Dec = boost::multiprecision::cpp_dec_float_100;

std::string format_sci_12(const Dec& x) {
    assert(x > 0);

    const Dec ten = 10;
    long long exp10 = floor(log10(x)).convert_to<long long>();
    Dec mant = x / pow(ten, exp10);

    const Dec scale = pow(ten, 12);
    mant = floor(mant * scale + Dec("0.5")) / scale;
    if (mant >= ten) {
        mant /= ten;
        ++exp10;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(12) << mant;
    oss << 'e' << exp10;
    return oss.str();
}

Dec E1_direct_m4() {
    const Dec q = Dec(1) - Dec(1) / 16;
    Dec sum = 0;

    for (int d = 1;; ++d) {
        const Dec qd = pow(q, d);
        const Dec term = Dec(d) * qd / (Dec(1) - qd);
        sum += term;
        if (term < Dec("1e-40")) {
            break;
        }
    }
    return sum;
}

Dec bernoulli_even(const int n) {
    if (n == 4 || n == 8) {
        return Dec(-1) / 30;
    }
    if (n == 16) {
        return Dec(-3617) / 510;
    }
    assert(false);
    return 0;
}

Dec E_odd_modular(const int k, const int p) {
    assert(k >= 3 && (k & 1) == 1);

    const int m = (k + 1) / 2;
    const int weight = 2 * m;

    const Dec q = Dec(1) - pow(Dec(2), -p);
    const Dec t = -log(q);
    const Dec y = t / (Dec(2) * boost::math::constants::pi<Dec>());

    const Dec B = bernoulli_even(weight);
    const Dec y_neg_pow = pow(y, -weight);

    // G = (1 - E_{2m}(iy)) * B_{2m} / (4m), and E_{2m}(iy) ~ y^{-2m}
    return (Dec(1) - y_neg_pow) * B / (Dec(4) * m);
}

}  // namespace

int main() {
    assert(format_sci_12(E1_direct_m4()) == "3.872155809243e2");
    assert(format_sci_12(E_odd_modular(3, 8)) == "2.767385314772e10");
    assert(format_sci_12(E_odd_modular(7, 15)) == "6.725803486744e39");

    std::cout << format_sci_12(E_odd_modular(15, 25)) << '\n';
    return 0;
}
