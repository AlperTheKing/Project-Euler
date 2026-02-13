#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace {

using boost::multiprecision::cpp_dec_float_100;

cpp_dec_float_100 f_value(int a, int b) {
    if (a == b) {
        return cpp_dec_float_100(0);
    }

    const int ab = a + b;
    cpp_dec_float_100 pow2 = 1;
    for (int i = 0; i < ab; ++i) {
        pow2 *= 2;
    }

    cpp_dec_float_100 term = 1;
    cpp_dec_float_100 sum = 1;

    for (int t = 0; t < 100000; ++t) {
        cpp_dec_float_100 ratio = 1;

        for (int i = 1; i <= ab; ++i) {
            ratio *= (ab * t + i);
        }
        for (int i = 1; i <= a; ++i) {
            ratio /= (a * t + i);
        }
        for (int i = 1; i <= b; ++i) {
            ratio /= (b * t + i);
        }
        ratio /= pow2;

        term *= ratio;
        sum += term;
        if (term < cpp_dec_float_100("1e-80")) {
            break;
        }
    }

    return cpp_dec_float_100(1) / sum;
}

void run_validations() {
    const cpp_dec_float_100 v11 = f_value(1, 1);
    assert(v11 == 0);

    const cpp_dec_float_100 v12 = f_value(1, 2);
    const cpp_dec_float_100 target = cpp_dec_float_100("0.427050983124842272");
    const cpp_dec_float_100 err = abs(v12 - target);
    assert(err < cpp_dec_float_100("1e-15"));
}

}  // namespace

int main() {
    run_validations();
    const cpp_dec_float_100 ans = f_value(89, 97);
    std::cout << std::fixed << std::setprecision(9) << ans << '\n';
    return 0;
}
