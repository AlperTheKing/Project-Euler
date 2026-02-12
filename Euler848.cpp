#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using Real = boost::multiprecision::cpp_dec_float_50;

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        s.push_back(static_cast<char>('0' + (x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u128 a_value(u64 m, u64 n) {
    u128 out = 0;
    if (n >= 1) out += std::min<u64>(m, 1);
    if (n >= 2) out += std::min<u64>(m, 2);
    if (n < 3) return out;

    out += std::min<u64>(m, 3);
    u64 idx = 4;
    u64 cnt = 3;
    u64 val = 6;

    while (idx <= n) {
        u64 take = std::min<u64>(cnt, n - idx + 1);
        out += static_cast<u128>(take) * std::min<u64>(m, val);
        idx += take;
        if (cnt > (std::numeric_limits<u64>::max() >> 1)) break;
        cnt <<= 1;
        if (val > (std::numeric_limits<u64>::max() >> 1)) break;
        val <<= 1;
    }
    return out;
}

static Real p_value(u64 m, u64 n) {
    u128 a = a_value(m, n);
    u128 den = static_cast<u128>(m) * n;
    Real ra(to_string_u128(a));
    Real rd(to_string_u128(den));
    return ra / rd;
}

int main() {
    Real sample = p_value(7, 5);
    assert(abs(sample - Real("0.51428571")) < Real("1e-8"));
    assert(p_value(1, 123) == Real(1));
    assert(abs(p_value(123, 1) - Real(1) / Real(123)) < Real("1e-30"));

    std::vector<u64> p7(21), p5(21);
    p7[0] = 1;
    p5[0] = 1;
    for (int i = 1; i <= 20; ++i) {
        p7[i] = p7[i - 1] * 7ULL;
        p5[i] = p5[i - 1] * 5ULL;
    }

    Real ans = 0;
    for (int i = 0; i <= 20; ++i) {
        for (int j = 0; j <= 20; ++j) {
            ans += p_value(p7[i], p5[j]);
        }
    }

    std::cout << std::fixed << std::setprecision(8) << ans << '\n';
    return 0;
}
