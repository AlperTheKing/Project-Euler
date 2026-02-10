#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using i64 = long long;

static std::vector<int> primes_upto(int n) {
    std::vector<bool> is_comp(n / 2 + 1, false);
    std::vector<int> primes;
    primes.push_back(2);
    for (int i = 3; (i64)i * i <= n; i += 2) {
        if (is_comp[i / 2]) continue;
        for (int j = i * i; j <= n; j += i * 2) is_comp[j / 2] = true;
    }
    for (int i = 3; i <= n; i += 2) {
        if (!is_comp[i / 2]) primes.push_back(i);
    }
    return primes;
}

static std::string fmt5(long double v) {
    assert(v > 0);
    int e = 0;
    long double m = v;
    while (m >= 10) {
        m /= 10;
        ++e;
    }
    while (m < 1) {
        m *= 10;
        --e;
    }

    const long double scale = 1e4L;
    long double mr = std::floor(m * scale + 0.5L) / scale;
    if (mr >= 10) {
        mr /= 10;
        ++e;
    }

    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(4) << mr << 'e' << e;
    return oss.str();
}

int main() {
    static constexpr int K = 7;
    static constexpr int P_MAX = 10'000'000;

    const auto primes = primes_upto(P_MAX);

    std::array<long double, K + 1> c{};
    c[0] = 1;
    for (int p : primes) {
        const long double a = 1.0L / ((long double)p * (long double)p);
        const long double b = 1 - a;
        for (int k = K; k >= 1; --k) c[k] = c[k] * b + c[k - 1] * a;
        c[0] *= b;
    }

    assert(fmt5(c[1]) == "3.3539e-1");
    assert(fmt5(c[2]) == "5.3293e-2");
    assert(fmt5(c[3]) == "3.2921e-3");
    assert(fmt5(c[4]) == "9.7046e-5");

    std::cout << fmt5(c[7]) << "\n";
    return 0;
}

