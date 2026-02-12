#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

using Poly = std::array<cpp_int, 6>;

static Poly mul_trunc(const Poly& a, const Poly& b) {
    Poly c{};
    for (int i = 0; i <= 5; ++i) {
        if (a[i] == 0) continue;
        for (int j = 0; j + i <= 5; ++j) {
            c[i + j] += a[i] * b[j];
        }
    }
    return c;
}

static Poly pow_poly_trunc(Poly base, std::uint64_t exp) {
    Poly res{};
    res[0] = 1;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = mul_trunc(res, base);
        }
        exp >>= 1ULL;
        if (exp > 0) {
            base = mul_trunc(base, base);
        }
    }
    return res;
}

static cpp_int pow_int(cpp_int base, std::uint64_t exp) {
    cpp_int res = 1;
    while (exp > 0) {
        if (exp & 1ULL) {
            res *= base;
        }
        base *= base;
        exp >>= 1ULL;
    }
    return res;
}

static cpp_int h_value(std::uint64_t m) {
    Poly R{};
    R[0] = 1;
    R[1] = 3;
    R[2] = 5;
    R[3] = 5;
    R[4] = 3;
    R[5] = 1;

    const Poly P = pow_poly_trunc(R, m);

    static const int B[6] = {1, 5, 10, 10, 5, 1};
    cpp_int h = 0;
    for (int i = 0; i <= 5; ++i) {
        h += P[i] * B[5 - i];
    }
    return h;
}

static std::string to_base7(cpp_int x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        cpp_int r = x % 7;
        s.push_back(static_cast<char>('0' + static_cast<int>(r)));
        x /= 7;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    const cpp_int h10 = h_value(10);
    const cpp_int g10 = h10 * pow_int(cpp_int(7), 10);
    assert(g10 == cpp_int("127278262644918"));

    const cpp_int h = h_value(142857);
    std::string s = to_base7(h);

    if (s.size() < 10) {
        s.append(10 - s.size(), '0');
    }

    std::cout << s.substr(0, 10) << '\n';
    return 0;
}
