#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

constexpr int kDigits = 10;

long double log10_factorial(int n) {
    static const long double kLn10 = std::logl(10.0L);
    return std::lgammal(static_cast<long double>(n) + 1.0L) / kLn10;
}

long double log10_s_value(int a, int b) {
    // S(a,b) = (b!)^a (a!)^b / (ab)!^2.
    long double log10_ab_fact = log10_factorial(a * b);
    return static_cast<long double>(a) * log10_factorial(b)
         + static_cast<long double>(b) * log10_factorial(a)
         - 2.0L * log10_ab_fact;
}

std::string format_scientific(long double log10_val, int digits) {
    long double exp_floor = std::floor(log10_val);
    long long exponent = static_cast<long long>(exp_floor);
    long double mantissa = std::powl(10.0L, log10_val - exp_floor);

    std::uint64_t scale = 1;
    for (int i = 0; i < digits; ++i) scale *= 10ULL;

    long double scaled_ld = mantissa * static_cast<long double>(scale);
    std::uint64_t scaled = static_cast<std::uint64_t>(std::floor(scaled_ld + 0.5L));

    if (scaled >= 10ULL * scale) {
        scaled /= 10ULL;
        exponent += 1;
    }

    std::uint64_t int_part = scaled / scale;
    std::uint64_t frac_part = scaled % scale;

    std::ostringstream oss;
    oss << int_part << '.' << std::setw(digits) << std::setfill('0') << frac_part
        << 'e' << exponent;
    return oss.str();
}

void validate() {
    struct TestCase {
        int a;
        int b;
        long double expected;
    };

    const TestCase tests[] = {
        {2, 2, 1.0L / 36.0L},
        {2, 3, 1.0L / 1800.0L},
    };

    for (const auto& test : tests) {
        long double log10_val = log10_s_value(test.a, test.b);
        long double val = std::powl(10.0L, log10_val);
        if (std::fabsl(val - test.expected) > 1e-15L) {
            std::cerr << "Validation failed for S(" << test.a << "," << test.b << ")." << '\n';
            std::exit(1);
        }
    }
}

}  // namespace

int main() {
    const int a = 5;
    const int b = 8;

    validate();

    long double log10_val = log10_s_value(a, b);
    std::string out = format_scientific(log10_val, kDigits);
    std::cout << out << '\n';
    return 0;
}
