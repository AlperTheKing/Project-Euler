#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetN = 10'000'000'000'000'000ULL;
constexpr u64 kMod = 282'475'249ULL; // 7^10

struct BeattySumsExact {
    u128 g = 0; // sum floor(k/phi)
    u128 p = 0; // sum k*floor(k/phi)
    u128 q = 0; // sum floor(k/phi)^2
};

struct BeattySumsMod {
    u64 g = 0;
    u64 p = 0;
    u64 q = 0;
};

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

u64 isqrt_u128(u128 n) {
    if (n == 0) {
        return 0;
    }

    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((static_cast<u128>(x) + 1) * (x + 1) <= n) {
        ++x;
    }
    while (static_cast<u128>(x) * x > n) {
        --x;
    }
    return x;
}

// floor(n / phi), computed exactly via floor(n*sqrt(5)).
u64 floor_div_phi(u64 n) {
    if (n == 0) {
        return 0;
    }

    const u128 nn = static_cast<u128>(n);
    const u128 radicand = static_cast<u128>(5) * nn * nn;
    const u64 floor_n_sqrt5 = isqrt_u128(radicand);
    return (floor_n_sqrt5 - n) / 2;
}

u128 triangular_u128(u64 n) {
    return static_cast<u128>(n) * (n + 1) / 2;
}

u128 square_sum_u128(u64 n) {
    return static_cast<u128>(n) * (n + 1) * (2 * static_cast<u128>(n) + 1) / 6;
}

u64 mod_add(u64 lhs, u64 rhs, u64 mod) {
    lhs += rhs;
    if (lhs >= mod) {
        lhs -= mod;
    }
    return lhs;
}

u64 mod_sub(u64 lhs, u64 rhs, u64 mod) {
    return (lhs >= rhs) ? (lhs - rhs) : (lhs + mod - rhs);
}

u64 mod_mul(u64 lhs, u64 rhs, u64 mod) {
    return static_cast<u64>((static_cast<u128>(lhs) * rhs) % mod);
}

u64 mod_inverse(u64 value, u64 mod) {
    i64 t = 0;
    i64 new_t = 1;
    i64 r = static_cast<i64>(mod);
    i64 new_r = static_cast<i64>(value % mod);

    while (new_r != 0) {
        const i64 quotient = r / new_r;

        const i64 next_t = t - quotient * new_t;
        t = new_t;
        new_t = next_t;

        const i64 next_r = r - quotient * new_r;
        r = new_r;
        new_r = next_r;
    }

    if (r != 1) {
        return 0;
    }

    if (t < 0) {
        t += static_cast<i64>(mod);
    }
    return static_cast<u64>(t);
}

u64 triangular_mod(u64 n, u64 mod, u64 inv2) {
    return mod_mul(mod_mul(n % mod, (n + 1) % mod, mod), inv2, mod);
}

u64 square_sum_mod(u64 n, u64 mod, u64 inv6) {
    const u64 a = n % mod;
    const u64 b = (n + 1) % mod;
    const u64 c = (2 * (n % mod) + 1) % mod;
    return mod_mul(mod_mul(mod_mul(a, b, mod), c, mod), inv6, mod);
}

BeattySumsExact beatty_sums_exact(u64 n) {
    if (n == 0) {
        return BeattySumsExact{};
    }

    const u64 m = floor_div_phi(n);
    const BeattySumsExact child = beatty_sums_exact(m);

    const u128 tri_m = triangular_u128(m);
    const u128 sq_m = square_sum_u128(m);
    const u128 n128 = static_cast<u128>(n);
    const u128 m128 = static_cast<u128>(m);

    const u128 g = n128 * m128 - tri_m - child.g;
    const u128 q = n128 * m128 * m128 - 2 * sq_m - 2 * child.p + tri_m + child.g;
    const u128 numerator = sq_m + 2 * child.p + child.q + tri_m + child.g;
    const u128 p = n128 * m128 * (n128 + 1) / 2 - numerator / 2;

    return BeattySumsExact{g, p, q};
}

BeattySumsMod beatty_sums_mod(u64 n, u64 mod, u64 inv2, u64 inv6) {
    if (n == 0) {
        return BeattySumsMod{};
    }

    const u64 m = floor_div_phi(n);
    const BeattySumsMod child = beatty_sums_mod(m, mod, inv2, inv6);

    const u64 n_mod = n % mod;
    const u64 m_mod = m % mod;
    const u64 tri_m = triangular_mod(m, mod, inv2);
    const u64 sq_m = square_sum_mod(m, mod, inv6);

    u64 g = mod_mul(n_mod, m_mod, mod);
    g = mod_sub(g, tri_m, mod);
    g = mod_sub(g, child.g, mod);

    u64 q = mod_mul(n_mod, mod_mul(m_mod, m_mod, mod), mod);
    q = mod_sub(q, mod_mul(2 % mod, sq_m, mod), mod);
    q = mod_sub(q, mod_mul(2 % mod, child.p, mod), mod);
    q = mod_add(q, tri_m, mod);
    q = mod_add(q, child.g, mod);

    u64 first = mod_mul(mod_mul(n_mod, m_mod, mod), (n + 1) % mod, mod);
    first = mod_mul(first, inv2, mod);

    u64 numerator = sq_m;
    numerator = mod_add(numerator, mod_mul(2 % mod, child.p, mod), mod);
    numerator = mod_add(numerator, child.q, mod);
    numerator = mod_add(numerator, tri_m, mod);
    numerator = mod_add(numerator, child.g, mod);

    const u64 second = mod_mul(numerator, inv2, mod);
    const u64 p = mod_sub(first, second, mod);

    return BeattySumsMod{g, p, q};
}

u128 solve_exact(u64 n) {
    const u64 cutoff = floor_div_phi(n + 1);
    const BeattySumsExact sums = beatty_sums_exact(cutoff);
    const u128 prefix = (4 * sums.p + sums.q + sums.g) / 2;

    if (cutoff == n) {
        return prefix;
    }

    const u128 n128 = static_cast<u128>(n);
    const u128 count = static_cast<u128>(n - cutoff);
    const u128 sum_x = triangular_u128(n) - triangular_u128(cutoff);
    const u128 sum_x2 = square_sum_u128(n) - square_sum_u128(cutoff);

    const u128 tail_numerator =
        count * (n128 * n128 + n128) + (2 * n128 - 1) * sum_x - 3 * sum_x2;
    const u128 tail = tail_numerator / 2;

    return prefix + tail;
}

u64 solve_mod(u64 n, u64 mod) {
    const u64 inv2 = mod_inverse(2, mod);
    const u64 inv6 = mod_inverse(6, mod);

    const u64 cutoff = floor_div_phi(n + 1);
    const BeattySumsMod sums = beatty_sums_mod(cutoff, mod, inv2, inv6);

    u64 prefix_numerator = mod_mul(4 % mod, sums.p, mod);
    prefix_numerator = mod_add(prefix_numerator, sums.q, mod);
    prefix_numerator = mod_add(prefix_numerator, sums.g, mod);
    const u64 prefix = mod_mul(prefix_numerator, inv2, mod);

    if (cutoff == n) {
        return prefix;
    }

    const u64 n_mod = n % mod;
    const u64 count_mod = (n - cutoff) % mod;
    const u64 sum_x = mod_sub(triangular_mod(n, mod, inv2),
                              triangular_mod(cutoff, mod, inv2),
                              mod);
    const u64 sum_x2 = mod_sub(square_sum_mod(n, mod, inv6),
                               square_sum_mod(cutoff, mod, inv6),
                               mod);

    const u64 n2_plus_n = mod_add(mod_mul(n_mod, n_mod, mod), n_mod, mod);
    const u64 coeff = mod_sub(mod_mul(2 % mod, n_mod, mod), 1 % mod, mod);

    const u64 term1 = mod_mul(count_mod, n2_plus_n, mod);
    const u64 term2 = mod_mul(coeff, sum_x, mod);
    const u64 term3 = mod_mul(3 % mod, sum_x2, mod);

    u64 tail_numerator = mod_add(term1, term2, mod);
    tail_numerator = mod_sub(tail_numerator, term3, mod);
    const u64 tail = mod_mul(tail_numerator, inv2, mod);

    return mod_add(prefix, tail, mod);
}

bool is_winning_configuration(u64 x, u64 y) {
    bool flip = false;

    while (true) {
        if (y % x == 0 || y / x >= 2) {
            return !flip;
        }

        y -= x; // only possible move when 1 < y/x < 2
        if (x > y) {
            std::swap(x, y);
        }
        flip = !flip;
    }
}

u128 brute_sum_by_game(u64 n, bool allow_multithreading, unsigned requested_threads = 0) {
    if (n < 2) {
        return 0;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    if (!allow_multithreading || n < 120 || threads <= 1) {
        threads = 1;
    } else {
        threads = std::min<unsigned>(threads, static_cast<unsigned>(n));
    }

    if (threads == 1) {
        u128 total = 0;
        for (u64 x = 1; x < n; ++x) {
            for (u64 y = x + 1; y <= n; ++y) {
                if (!is_winning_configuration(x, y)) {
                    total += x + y;
                }
            }
        }
        return total;
    }

    std::atomic<u64> next_x{1};
    std::vector<u128> partial(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            u128 local = 0;
            while (true) {
                const u64 x = next_x.fetch_add(1, std::memory_order_relaxed);
                if (x >= n) {
                    break;
                }
                for (u64 y = x + 1; y <= n; ++y) {
                    if (!is_winning_configuration(x, y)) {
                        local += x + y;
                    }
                }
            }
            partial[tid] = local;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u128 total = 0;
    for (const u128 value : partial) {
        total += value;
    }
    return total;
}

bool run_validation_checkpoints(bool allow_multithreading, unsigned requested_threads) {
    {
        constexpr u64 kBruteN = 220;
        const u128 brute = brute_sum_by_game(kBruteN, allow_multithreading, requested_threads);
        const u128 fast = solve_exact(kBruteN);
        if (brute != fast) {
            std::cerr << "Brute-force checkpoint failed for N=" << kBruteN << ": got "
                      << to_string_u128(fast) << ", expected " << to_string_u128(brute)
                      << "\n";
            return false;
        }
    }

    {
        constexpr u64 kN = 10;
        constexpr u64 kExpected = 211;
        const u128 got = solve_exact(kN);
        if (got != kExpected) {
            std::cerr << "Checkpoint failed for S(10): got " << to_string_u128(got)
                      << ", expected " << kExpected << "\n";
            return false;
        }
    }

    {
        constexpr u64 kN = 10'000;
        constexpr u64 kExpected = 230'312'207'313ULL;
        const u128 got = solve_exact(kN);
        if (got != kExpected) {
            std::cerr << "Checkpoint failed for S(10^4): got " << to_string_u128(got)
                      << ", expected " << kExpected << "\n";
            return false;
        }
    }

    {
        constexpr u64 kN = 1'000'000;
        const u128 exact = solve_exact(kN);
        const u64 exact_mod = static_cast<u64>(exact % kMod);
        const u64 fast_mod = solve_mod(kN, kMod);
        if (fast_mod != exact_mod) {
            std::cerr << "Modulo checkpoint failed for N=" << kN << ": got " << fast_mod
                      << ", expected " << exact_mod << "\n";
            return false;
        }
    }

    return true;
}

} // namespace

int main() {
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    if (!run_validation_checkpoints(true, threads)) {
        return 1;
    }

    const u64 answer = solve_mod(kTargetN, kMod);
    std::cout << answer << '\n';
    return 0;
}
