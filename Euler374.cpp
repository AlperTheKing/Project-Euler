#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr i64 kMod = 982451653LL;
constexpr u64 kLimit = 100000000000000ULL;

u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

u64 max_m_for_n(const u64 n) {
    // Largest m with m(m+3)/2 <= n.
    u64 m = (isqrt_u64(8ULL * n + 9ULL) - 3ULL) / 2ULL;
    while ((m + 1ULL) * (m + 4ULL) / 2ULL <= n) {
        ++m;
    }
    while (m > 0ULL && m * (m + 3ULL) / 2ULL > n) {
        --m;
    }
    return m;
}

u128 factorial_u128(const u64 n) {
    u128 f = 1U;
    for (u64 i = 2; i <= n; ++i) {
        f *= static_cast<u128>(i);
    }
    return f;
}

u128 value_exact_small(const u64 n) {
    if (n == 0ULL) {
        return 0U;
    }
    if (n == 1ULL) {
        return 1U;
    }

    const u64 m = max_m_for_n(n);
    const u64 t = m * (m + 3ULL) / 2ULL;
    const u64 s = n - t;

    const u128 fact_m2 = factorial_u128(m + 2ULL);
    if (s <= m - 1ULL) {
        return static_cast<u128>(m) * (fact_m2 / static_cast<u128>(m + 2ULL - s));
    }
    if (s == m) {
        return static_cast<u128>(m) * (fact_m2 / 2U);
    }
    return static_cast<u128>(m) * (factorial_u128(m + 3ULL) / (2U * static_cast<u128>(m + 2ULL)));
}

i64 sum_mod(const u64 n) {
    if (n == 0ULL) {
        return 0;
    }

    i64 answer = 1;  // n = 1
    const u64 m_max = max_m_for_n(n);
    const u64 inv_limit = m_max + 3ULL;

    std::vector<i64> inv(static_cast<std::size_t>(inv_limit + 1ULL), 0);
    inv[1] = 1;
    for (u64 i = 2ULL; i <= inv_limit; ++i) {
        inv[static_cast<std::size_t>(i)] =
            (kMod - static_cast<i64>((static_cast<__int128>(kMod / static_cast<i64>(i)) *
                                      inv[static_cast<std::size_t>(kMod % static_cast<i64>(i))]) %
                                     kMod)) %
            kMod;
    }

    const i64 inv2 = (kMod + 1LL) / 2LL;

    i64 fact = 1;          // running factorial up to `upto`
    i64 harmonic = 0;      // sum_{d=2..upto} inv[d] mod kMod
    u64 upto = 1ULL;

    for (u64 m = 1ULL; m <= m_max; ++m) {
        const u64 target = m + 2ULL;
        while (upto < target) {
            ++upto;
            fact = static_cast<i64>((static_cast<__int128>(fact) * static_cast<i64>(upto)) % kMod);
            if (upto >= 2ULL) {
                harmonic += inv[static_cast<std::size_t>(upto)];
                if (harmonic >= kMod) {
                    harmonic -= kMod;
                }
            }
        }

        const u64 start = m * (m + 3ULL) / 2ULL;
        const u64 end_full = (m + 1ULL) * (m + 4ULL) / 2ULL - 1ULL;
        const u64 len = std::min(n, end_full) - start + 1ULL;

        const i64 mm = static_cast<i64>(m % static_cast<u64>(kMod));
        const i64 base = static_cast<i64>((static_cast<__int128>(mm) * fact) % kMod);

        i64 add = 0;
        if (len <= m) {
            const u64 d1 = m + 3ULL - len;
            const u64 d2 = m + 2ULL;
            i64 partial_h = 0;
            for (u64 d = d1; d <= d2; ++d) {
                partial_h += inv[static_cast<std::size_t>(d)];
                if (partial_h >= kMod) {
                    partial_h -= kMod;
                }
            }
            add = static_cast<i64>((static_cast<__int128>(base) * partial_h) % kMod);
        } else if (len == m + 1ULL) {
            add = static_cast<i64>((static_cast<__int128>(base) * harmonic) % kMod);
        } else {  // len == m + 2
            add = static_cast<i64>((static_cast<__int128>(base) * harmonic) % kMod);
            const i64 fact_m3 = static_cast<i64>((static_cast<__int128>(fact) * static_cast<i64>(m + 3ULL)) % kMod);
            i64 extra = static_cast<i64>((static_cast<__int128>(mm) * fact_m3) % kMod);
            extra = static_cast<i64>((static_cast<__int128>(extra) * inv2) % kMod);
            extra = static_cast<i64>((static_cast<__int128>(extra) * inv[static_cast<std::size_t>(m + 2ULL)]) % kMod);
            add += extra;
            if (add >= kMod) {
                add -= kMod;
            }
        }

        answer += add;
        if (answer >= kMod) {
            answer -= kMod;
        }
    }

    return answer;
}

bool run_checkpoints() {
    if (value_exact_small(5ULL) != 12U) {
        std::cerr << "Checkpoint failed: f(5)*m(5)\n";
        return false;
    }
    if (value_exact_small(10ULL) != 90U) {
        std::cerr << "Checkpoint failed: f(10)*m(10)\n";
        return false;
    }

    u128 sum_100 = 0U;
    for (u64 n = 1ULL; n <= 100ULL; ++n) {
        sum_100 += value_exact_small(n);
    }
    if (sum_100 != 1683550844462ULL) {
        std::cerr << "Checkpoint failed: sum up to 100\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::cout << sum_mod(kLimit) << '\n';
    return 0;
}
