#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;

constexpr u64 TARGET_M = 1234ULL * 100'000'000ULL;
constexpr u64 TARGET_N = 2345ULL * 100'000'000ULL;
constexpr u32 MOD = 1'000'000'007U;
constexpr int SIEVE_LIMIT = 5'000'000;
constexpr u64 INV2 = 500'000'004ULL;
constexpr u64 INV3 = 333'333'336ULL;
constexpr u64 INV6 = 166'666'668ULL;

u32 add_mod(const u32 a, const u32 b) {
    const u32 s = a + b;
    return s >= MOD ? s - MOD : s;
}

u32 sub_mod(const u32 a, const u32 b) {
    return a >= b ? a - b : static_cast<u32>(static_cast<u64>(a) + MOD - b);
}

u32 mul_mod(const u64 a, const u64 b) {
    return static_cast<u32>((a % MOD) * (b % MOD) % MOD);
}

u32 sum_power(const int power, const u64 n) {
    const u64 a = n % MOD;
    const u64 b = (n + 1) % MOD;
    if (power == 0) {
        return static_cast<u32>(a);
    }
    if (power == 1) {
        return static_cast<u32>(a * b % MOD * INV2 % MOD);
    }
    if (power == 2) {
        const u64 c = (2 * (n % MOD) + 1) % MOD;
        return static_cast<u32>(a * b % MOD * c % MOD * INV6 % MOD);
    }
    const u64 t = a * b % MOD * INV2 % MOD;
    return static_cast<u32>(t * t % MOD);
}

u32 range_power(const int power, const u64 lo, const u64 hi) {
    return sub_mod(sum_power(power, hi), sum_power(power, lo - 1));
}

u32 choose2(const u64 n) {
    const u64 a = n % MOD;
    const u64 b = (n - 1) % MOD;
    return static_cast<u32>(a * b % MOD * INV2 % MOD);
}

u32 choose3(const u64 n) {
    if (n < 3) {
        return 0;
    }
    const u64 a = n % MOD;
    const u64 b = (n - 1) % MOD;
    const u64 c = (n - 2) % MOD;
    return static_cast<u32>(a * b % MOD * c % MOD * INV6 % MOD);
}

u32 triangle_number(const u64 n) {
    const u64 a = n % MOD;
    const u64 b = (n + 1) % MOD;
    return static_cast<u32>(a * b % MOD * INV2 % MOD);
}

struct TotientSummatory {
    explicit TotientSummatory(const int limit) : limit(limit) {
        build();
        for (auto& table : memo) {
            table.reserve(200'000);
        }
    }

    u32 prefix(const int power, const u64 n) {
        if (n <= static_cast<u64>(limit)) {
            return sums[static_cast<std::size_t>(power)][static_cast<std::size_t>(n)];
        }

        auto& table = memo[static_cast<std::size_t>(power)];
        const auto it = table.find(n);
        if (it != table.end()) {
            return it->second;
        }

        u32 result = sum_power(power + 1, n);
        for (u64 lo = 2; lo <= n;) {
            const u64 q = n / lo;
            const u64 hi = n / q;
            const u32 block = range_power(power, lo, hi);
            result = sub_mod(result, mul_mod(block, prefix(power, q)));
            lo = hi + 1;
        }

        table.emplace(n, result);
        return result;
    }

private:
    int limit;
    std::array<std::vector<u32>, 3> sums;
    std::array<std::unordered_map<u64, u32>, 3> memo;

    void build() {
        std::vector<int> phi(static_cast<std::size_t>(limit + 1), 0);
        std::vector<int> primes;
        std::vector<unsigned char> composite(static_cast<std::size_t>(limit + 1), 0);
        phi[1] = 1;

        for (int i = 2; i <= limit; ++i) {
            if (!composite[static_cast<std::size_t>(i)]) {
                primes.push_back(i);
                phi[static_cast<std::size_t>(i)] = i - 1;
            }
            for (const int p : primes) {
                const i64 v = static_cast<i64>(i) * p;
                if (v > limit) {
                    break;
                }
                composite[static_cast<std::size_t>(v)] = 1;
                if (i % p == 0) {
                    phi[static_cast<std::size_t>(v)] = phi[static_cast<std::size_t>(i)] * p;
                    break;
                }
                phi[static_cast<std::size_t>(v)] = phi[static_cast<std::size_t>(i)] * (p - 1);
            }
        }

        for (auto& row : sums) {
            row.assign(static_cast<std::size_t>(limit + 1), 0);
        }
        for (int i = 1; i <= limit; ++i) {
            const u64 im = static_cast<u64>(i) % MOD;
            const u64 ph = static_cast<u64>(phi[static_cast<std::size_t>(i)]) % MOD;
            sums[0][static_cast<std::size_t>(i)] = add_mod(sums[0][static_cast<std::size_t>(i - 1)],
                                                           static_cast<u32>(ph));
            sums[1][static_cast<std::size_t>(i)] = add_mod(sums[1][static_cast<std::size_t>(i - 1)],
                                                           static_cast<u32>(im * ph % MOD));
            sums[2][static_cast<std::size_t>(i)] = add_mod(sums[2][static_cast<std::size_t>(i - 1)],
                                                           static_cast<u32>(im * im % MOD * ph % MOD));
        }
    }
};

u32 pairwise_triangle_candidates(const u64 m, const u64 n) {
    u32 total = mul_mod(choose3(m), choose3(n));
    total = add_mod(total, mul_mod(mul_mod(choose3(m), n % MOD), (n - 1) % MOD));
    total = add_mod(total, mul_mod(mul_mod(m % MOD, (m - 1) % MOD), choose3(n)));
    const u32 corner = mul_mod(mul_mod(m % MOD, (m - 1) % MOD), mul_mod(n % MOD, (n - 1) % MOD));
    total = add_mod(total, mul_mod(corner, INV2));
    return total;
}

u32 concurrent_triples(const u64 m, const u64 n, TotientSummatory& totient) {
    const u64 limit = std::min(m, n) - 1;
    u32 total = 0;

    for (u64 lo = 2; lo <= limit;) {
        const u64 qm = (m - 1) / lo;
        const u64 qn = (n - 1) / lo;
        const u64 hi = std::min((m - 1) / qm, (n - 1) / qn);

        const u32 s0 = sub_mod(totient.prefix(0, hi), totient.prefix(0, lo - 1));
        const u32 s1 = sub_mod(totient.prefix(1, hi), totient.prefix(1, lo - 1));
        const u32 s2 = sub_mod(totient.prefix(2, hi), totient.prefix(2, lo - 1));

        const u32 am0 = mul_mod(qm % MOD, m % MOD);
        const u32 an0 = mul_mod(qn % MOD, n % MOD);
        const u32 am1 = triangle_number(qm);
        const u32 an1 = triangle_number(qn);

        u32 term = mul_mod(mul_mod(am0, an0), s0);
        term = sub_mod(term, mul_mod(add_mod(mul_mod(am0, an1), mul_mod(am1, an0)), s1));
        term = add_mod(term, mul_mod(mul_mod(am1, an1), s2));
        total = add_mod(total, term);

        lo = hi + 1;
    }

    return total;
}

u32 solve_mod(const u64 m, const u64 n, TotientSummatory& totient) {
    return sub_mod(pairwise_triangle_candidates(m, n), concurrent_triples(m, n, totient));
}

bool pair_intersects(const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return a.first == b.first || a.second == b.second ||
           static_cast<i64>(a.first - b.first) * static_cast<i64>(a.second - b.second) < 0;
}

bool collinear(const std::pair<int, int>& a, const std::pair<int, int>& b, const std::pair<int, int>& c) {
    return static_cast<i64>(b.first - a.first) * static_cast<i64>(c.second - a.second) ==
           static_cast<i64>(b.second - a.second) * static_cast<i64>(c.first - a.first);
}

u64 brute_count(const int m, const int n) {
    std::vector<std::pair<int, int>> segments;
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            segments.emplace_back(i, j);
        }
    }

    u64 total = 0;
    for (std::size_t a = 0; a < segments.size(); ++a) {
        for (std::size_t b = a + 1; b < segments.size(); ++b) {
            for (std::size_t c = b + 1; c < segments.size(); ++c) {
                if (pair_intersects(segments[a], segments[b]) &&
                    pair_intersects(segments[a], segments[c]) &&
                    pair_intersects(segments[b], segments[c]) &&
                    !collinear(segments[a], segments[b], segments[c])) {
                    ++total;
                }
            }
        }
    }
    return total;
}

void run_checkpoints(TotientSummatory& totient) {
    for (int m = 1; m <= 5; ++m) {
        for (int n = 1; n <= 5; ++n) {
            assert(solve_mod(m, n, totient) == brute_count(m, n));
        }
    }

    assert(solve_mod(2, 3, totient) == 8U);
    assert(solve_mod(3, 5, totient) == 146U);
    assert(solve_mod(12, 23, totient) == 756'716U);
}

}  // namespace

int main(int argc, char** argv) {
    bool should_run_checkpoints = true;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--skip-checkpoints") {
            should_run_checkpoints = false;
        }
    }

    TotientSummatory totient(SIEVE_LIMIT);
    if (should_run_checkpoints) {
        run_checkpoints(totient);
    }

    std::cout << solve_mod(TARGET_M, TARGET_N, totient) << '\n';
    return 0;
}
