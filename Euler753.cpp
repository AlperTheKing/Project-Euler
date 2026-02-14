#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

i64 isqrt_i64(const i64 n) {
    long double d = static_cast<long double>(n);
    i64 x = static_cast<i64>(std::sqrt(d));
    while ((x + 1) * (x + 1) <= n) {
        ++x;
    }
    while (x * x > n) {
        --x;
    }
    return x;
}

i64 F_formula(const int p, const std::vector<i64>& twenty_seven_u2) {
    if (p == 3) {
        return 2;
    }
    if (p < 5 || p % 3 == 0) {
        return 0;
    }
    if (p % 3 == 2) {
        return static_cast<i64>(p - 1) * static_cast<i64>(p - 2);
    }

    const i64 target = 4LL * p;
    i64 t_sel = 0;
    for (const i64 v : twenty_seven_u2) {
        if (v > target) {
            break;
        }
        const i64 rem = target - v;
        const i64 t = isqrt_i64(rem);
        if (t * t == rem) {
            t_sel = (t % 3 == 1) ? t : -t;
            break;
        }
    }
    assert(t_sel != 0);

    return static_cast<i64>(p - 1) * static_cast<i64>(p - 8 + t_sel);
}

i64 F_bruteforce(const int p) {
    i64 count = 0;
    std::vector<int> cubes(static_cast<std::size_t>(p), 0);
    for (int x = 0; x < p; ++x) {
        cubes[static_cast<std::size_t>(x)] = static_cast<int>((1LL * x * x * x) % p);
    }

    for (int a = 1; a < p; ++a) {
        const int a3 = cubes[static_cast<std::size_t>(a)];
        for (int b = 1; b < p; ++b) {
            const int s = (a3 + cubes[static_cast<std::size_t>(b)]) % p;
            for (int c = 1; c < p; ++c) {
                if (cubes[static_cast<std::size_t>(c)] == s) {
                    ++count;
                }
            }
        }
    }
    return count;
}

i64 solve() {
    constexpr int kLimit = 6'000'000;
    std::vector<bool> is_prime(static_cast<std::size_t>(kLimit), true);
    is_prime[0] = false;
    is_prime[1] = false;
    for (int i = 2; 1LL * i * i < kLimit; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j < kLimit; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<i64> twenty_seven_u2;
    for (i64 u = 1;; ++u) {
        const i64 v = 27LL * u * u;
        if (v > 4LL * kLimit) {
            break;
        }
        twenty_seven_u2.push_back(v);
    }

    i64 total = 0;
    for (int p = 2; p < kLimit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        total += F_formula(p, twenty_seven_u2);
    }

    return total;
}

}  // namespace

int main() {
    std::vector<i64> test_u2;
    for (i64 u = 1; 27LL * u * u <= 4LL * 10'000; ++u) {
        test_u2.push_back(27LL * u * u);
    }

    assert(F_formula(5, test_u2) == 12);
    assert(F_formula(7, test_u2) == 0);

    for (int p = 5; p < 200; ++p) {
        bool prime = true;
        for (int d = 2; d * d <= p; ++d) {
            if (p % d == 0) {
                prime = false;
                break;
            }
        }
        if (!prime) {
            continue;
        }
        assert(F_formula(p, test_u2) == F_bruteforce(p));
    }

    std::cout << solve() << '\n';
    return 0;
}
