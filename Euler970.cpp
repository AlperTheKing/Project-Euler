#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Dec = boost::multiprecision::cpp_dec_float_100;

std::string extract_non_six_digits_from_fraction(Dec fractional_value, int wanted) {
    std::string out;
    out.reserve(static_cast<std::size_t>(wanted));

    Dec f = fractional_value;
    if (f < 0) {
        f = f - boost::multiprecision::floor(f);
    }

    for (int step = 0; static_cast<int>(out.size()) < wanted && step < 100000; ++step) {
        f *= 10;
        int digit = static_cast<int>(boost::multiprecision::floor(f + Dec("1e-70")));
        if (digit < 0) {
            digit = 0;
        }
        if (digit > 9) {
            digit = 9;
        }
        f -= Dec(digit);
        if (digit != 6) {
            out.push_back(static_cast<char>('0' + digit));
        }
    }

    return out;
}

std::string non_six_digits_large_n(std::int64_t n, int wanted) {
    // Dominant nonzero pole of 1 / (s - 1 + e^{-s}), i.e. 1 + W_1(-1/e).
    const Dec alpha("-2.0888430156130438559570867167749475005456937410367296732391125442446071101031945");
    const Dec beta("7.4614892856542545569061166121864153345090949932022092409344113914118766543223747");
    const Dec den = alpha * alpha + beta * beta;

    const Dec two_pi = 2 * boost::math::constants::pi<Dec>();
    const Dec raw = beta * Dec(n);
    const Dec theta = raw - boost::multiprecision::floor(raw / two_pi) * two_pi;
    const Dec c = boost::multiprecision::cos(theta);
    const Dec s = boost::multiprecision::sin(theta);

    const Dec envelope = 2 * (c * alpha + s * beta) / den;
    const Dec sign = (envelope >= 0 ? Dec(1) : Dec(-1));

    const Dec log10_eps =
        alpha * Dec(n) / boost::multiprecision::log(Dec(10)) +
        boost::multiprecision::log10(boost::multiprecision::abs(envelope));

    const std::int64_t shift = static_cast<std::int64_t>(
        boost::multiprecision::floor(-log10_eps).convert_to<long long>());
    const Dec frac = -log10_eps - Dec(shift);
    const Dec t = sign * boost::multiprecision::pow(Dec(10), frac);  // epsilon * 10^shift

    const Dec base = Dec(2) / 3 + t;
    const int q = static_cast<int>(boost::multiprecision::floor(base));
    Dec r = base - q;
    if (r < 0) {
        r += 1;
    }

    // Last few digits before the shifted position: ...666666 + q.
    std::vector<int> tail(40, 6);
    int carry = q;
    for (int i = static_cast<int>(tail.size()) - 1; i >= 0 && carry != 0; --i) {
        const int value = tail[static_cast<std::size_t>(i)] + carry;
        if (value >= 0) {
            tail[static_cast<std::size_t>(i)] = value % 10;
            carry = value / 10;
        } else {
            const int borrow = (-value + 9) / 10;
            tail[static_cast<std::size_t>(i)] = value + 10 * borrow;
            carry = -borrow;
        }
    }

    std::string out;
    out.reserve(static_cast<std::size_t>(wanted));

    for (int d : tail) {
        if (d != 6) {
            out.push_back(static_cast<char>('0' + d));
            if (static_cast<int>(out.size()) == wanted) {
                return out;
            }
        }
    }

    std::string rest = extract_non_six_digits_from_fraction(r, wanted - static_cast<int>(out.size()));
    out += rest;
    return out;
}

bool run_checkpoints() {
    const Dec e = boost::multiprecision::exp(Dec(1));
    const Dec h2 = e * e - e;
    const Dec h3 = e * e * e - 2 * e * e + e / 2;

    const std::string d2 =
        extract_non_six_digits_from_fraction(h2 - boost::multiprecision::floor(h2), 8);
    const std::string d3 =
        extract_non_six_digits_from_fraction(h3 - boost::multiprecision::floor(h3), 8);

    if (d2 != "70774270") {
        std::cerr << "Checkpoint failed for H(2): got " << d2 << '\n';
        return false;
    }
    if (d3 != "55395558") {
        std::cerr << "Checkpoint failed for H(3): got " << d3 << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    std::cout << non_six_digits_large_n(1'000'000, 8) << '\n';
    return 0;
}
