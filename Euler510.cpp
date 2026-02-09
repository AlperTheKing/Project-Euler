#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

u64 isqrt_u128(const u128 x) {
    long double approx = std::sqrt(static_cast<long double>(x));
    u64 r = static_cast<u64>(approx);
    while (static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > x) {
        --r;
    }
    return r;
}

u64 sum_formula(const u64 n_max) {
    u128 total = 0;

    for (u64 n = 1ULL;; ++n) {
        const u128 min_b = static_cast<u128>(n + 1ULL) * static_cast<u128>(n + 1ULL) *
                           static_cast<u128>(n) * static_cast<u128>(n);
        if (min_b > n_max) {
            break;
        }

        for (u64 m = 1ULL; m <= n; ++m) {
            if (std::gcd(m, n) != 1ULL) {
                continue;
            }

            const u128 sum = static_cast<u128>(m + n);
            const u128 b0 = sum * sum * static_cast<u128>(n) * static_cast<u128>(n);
            if (b0 > n_max) {
                break;
            }
            const u128 a0 = sum * sum * static_cast<u128>(m) * static_cast<u128>(m);
            const u128 c0 = static_cast<u128>(m) * static_cast<u128>(m) * static_cast<u128>(n) *
                            static_cast<u128>(n);

            const u128 t_max = static_cast<u128>(n_max) / b0;
            const u128 t_sum = t_max * (t_max + 1ULL) / 2ULL;
            total += (a0 + b0 + c0) * t_sum;
        }
    }

    return static_cast<u64>(total);
}

u64 sum_bruteforce(const u64 n_max) {
    u128 total = 0;
    for (u64 a = 1ULL; a <= n_max; ++a) {
        for (u64 b = a; b <= n_max; ++b) {
            const u128 prod = static_cast<u128>(a) * static_cast<u128>(b);
            const u64 s = isqrt_u128(prod);
            if (static_cast<u128>(s) * static_cast<u128>(s) != prod) {
                continue;
            }

            const u128 denom = static_cast<u128>(a + b) + 2ULL * static_cast<u128>(s);
            if (prod % denom != 0ULL) {
                continue;
            }

            const u128 c = prod / denom;
            total += static_cast<u128>(a + b) + c;
        }
    }
    return static_cast<u64>(total);
}

bool run_checkpoints() {
    if (sum_formula(5ULL) != 9ULL) {
        std::cerr << "Checkpoint failed: S(5)\n";
        return false;
    }
    if (sum_formula(100ULL) != 3'072ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    if (sum_formula(300ULL) != sum_bruteforce(300ULL)) {
        std::cerr << "Checkpoint failed: formula/bruteforce mismatch at 300\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 1'000'000'000ULL;
    std::cout << sum_formula(n) << '\n';
    return 0;
}
