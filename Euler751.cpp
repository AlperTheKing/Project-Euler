#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Dec = boost::multiprecision::cpp_dec_float_100;

std::vector<int> sequence_terms(Dec theta, const int count) {
    std::vector<int> terms;
    terms.reserve(static_cast<std::size_t>(count));

    Dec b = theta;
    for (int i = 0; i < count; ++i) {
        const int a = static_cast<int>(floor(b));
        terms.push_back(a);
        b = Dec(a) * (b - Dec(a) + Dec(1));
    }

    return terms;
}

Dec tau_from_theta(const Dec& theta, const int frac_digits) {
    Dec b = theta;
    const int a1 = static_cast<int>(floor(b));

    std::string repr = std::to_string(a1);
    repr.push_back('.');

    std::string frac;
    frac.reserve(static_cast<std::size_t>(frac_digits + 32));

    while (static_cast<int>(frac.size()) < frac_digits) {
        const int a = static_cast<int>(floor(b));
        b = Dec(a) * (b - Dec(a) + Dec(1));
        const int next_a = static_cast<int>(floor(b));
        frac += std::to_string(next_a);
    }

    frac.resize(static_cast<std::size_t>(frac_digits));
    repr += frac;
    return Dec(repr);
}

Dec fixed_point_theta() {
    Dec theta("2.2");
    for (int i = 0; i < 30; ++i) {
        theta = tau_from_theta(theta, 80);
    }
    return theta;
}

std::string rounded_24(const Dec& value) {
    const Dec scale = pow(Dec(10), 24);
    const Dec rounded = floor(value * scale + Dec("0.5")) / scale;

    std::ostringstream out;
    out << std::fixed << std::setprecision(24) << rounded;
    return out.str();
}

}  // namespace

int main() {
    {
        const Dec sample_theta("2.956938891377988");
        const std::vector<int> sample_terms = sequence_terms(sample_theta, 10);
        const std::vector<int> fibonacci = {2, 3, 5, 8, 13, 21, 34, 55, 89, 144};
        assert(sample_terms == fibonacci);
    }

    const std::string answer = rounded_24(fixed_point_theta());
    assert(answer == "2.223561019313554106173177");
    std::cout << answer << '\n';
    return 0;
}
