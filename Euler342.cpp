#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u64 integer_cuberoot_floor(const u128 x) {
    long double approx = std::cbrt(static_cast<long double>(x));
    u64 r = static_cast<u64>(approx);

    while (static_cast<u128>(r + 1) * (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * r * r > x) {
        --r;
    }
    return r;
}

bool is_perfect_cube(const u128 x) {
    const u64 r = integer_cuberoot_floor(x);
    return static_cast<u128>(r) * r * r == x;
}

u128 brute_sum(const int limit) {
    std::vector<int> phi(static_cast<std::size_t>(limit));
    for (int i = 0; i < limit; ++i) {
        phi[static_cast<std::size_t>(i)] = i;
    }

    for (int p = 2; p < limit; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (int x = p; x < limit; x += p) {
            phi[static_cast<std::size_t>(x)] -= phi[static_cast<std::size_t>(x)] / p;
        }
    }

    u128 sum = 0;
    for (int n = 2; n < limit; ++n) {
        const u128 value = static_cast<u128>(n) * static_cast<u128>(phi[static_cast<std::size_t>(n)]);
        if (is_perfect_cube(value)) {
            sum += static_cast<u128>(n);
        }
    }
    return sum;
}

struct Solver {
    u64 limit_exclusive;
    int max_prime;

    std::vector<int> primes_desc;
    std::vector<int> spf;
    std::vector<std::vector<std::pair<int, int>>> factors_mod3_of_p_minus_1;

    std::vector<std::uint8_t> residue_mod3;
    std::vector<int> included_primes;

    u128 answer = 0;

    explicit Solver(const u64 limit)
        : limit_exclusive(limit), max_prime(static_cast<int>(std::sqrt(static_cast<long double>(limit - 1)))) {}

    void build_primes_and_spf() {
        spf.assign(static_cast<std::size_t>(max_prime + 1), 0);
        for (int i = 0; i <= max_prime; ++i) {
            spf[static_cast<std::size_t>(i)] = i;
        }

        for (int i = 2; static_cast<i64>(i) * i <= max_prime; ++i) {
            if (spf[static_cast<std::size_t>(i)] != i) {
                continue;
            }
            for (int j = i * i; j <= max_prime; j += i) {
                if (spf[static_cast<std::size_t>(j)] == j) {
                    spf[static_cast<std::size_t>(j)] = i;
                }
            }
        }

        std::vector<int> primes_asc;
        for (int i = 2; i <= max_prime; ++i) {
            if (spf[static_cast<std::size_t>(i)] == i) {
                primes_asc.push_back(i);
            }
        }

        primes_desc = primes_asc;
        std::reverse(primes_desc.begin(), primes_desc.end());

        factors_mod3_of_p_minus_1.assign(static_cast<std::size_t>(max_prime + 1), {});
        for (int p : primes_asc) {
            int x = p - 1;
            auto& vec = factors_mod3_of_p_minus_1[static_cast<std::size_t>(p)];
            while (x > 1) {
                const int q = spf[static_cast<std::size_t>(x)];
                int cnt = 0;
                while (x % q == 0) {
                    x /= q;
                    ++cnt;
                }
                cnt %= 3;
                if (cnt != 0) {
                    vec.emplace_back(q, cnt);
                }
            }
        }

        residue_mod3.assign(static_cast<std::size_t>(max_prime + 1), 0);
    }

    void apply_factors(const int p, const int sign) {
        for (const auto& [q, e] : factors_mod3_of_p_minus_1[static_cast<std::size_t>(p)]) {
            int v = static_cast<int>(residue_mod3[static_cast<std::size_t>(q)]);
            if (sign > 0) {
                v += e;
            } else {
                v -= e;
            }
            v %= 3;
            if (v < 0) {
                v += 3;
            }
            residue_mod3[static_cast<std::size_t>(q)] = static_cast<std::uint8_t>(v);
        }
    }

    u64 pow_u64(const u64 base, const int exp) const {
        u64 value = 1;
        for (int i = 0; i < exp; ++i) {
            value *= base;
        }
        return value;
    }

    u128 multiplier_sum_rec(const int pos, const u64 remaining) const {
        if (pos == static_cast<int>(included_primes.size())) {
            return 1;
        }

        const u64 p = static_cast<u64>(included_primes[static_cast<std::size_t>(pos)]);
        const u128 p3 = static_cast<u128>(p) * p * p;

        u128 total = 0;
        u128 mul = 1;

        while (mul <= static_cast<u128>(remaining)) {
            const u64 next_remaining = remaining / static_cast<u64>(mul);
            total += mul * multiplier_sum_rec(pos + 1, next_remaining);

            if (mul > static_cast<u128>(remaining) / p3) {
                break;
            }
            mul *= p3;
        }

        return total;
    }

    void add_solution(const u64 base_n) {
        if (base_n <= 1 || base_n >= limit_exclusive) {
            return;
        }

        const u64 remaining = (limit_exclusive - 1) / base_n;
        const u128 mul_sum = multiplier_sum_rec(0, remaining);
        answer += static_cast<u128>(base_n) * mul_sum;
    }

    void dfs(int idx, const u64 current_n) {
        int i = idx;
        const int npr = static_cast<int>(primes_desc.size());

        // Skip deterministic exclusions: residue is zero and p^2 already too large.
        while (i < npr) {
            const int p = primes_desc[static_cast<std::size_t>(i)];
            const int c = residue_mod3[static_cast<std::size_t>(p)];
            if (c != 0) {
                break;
            }

            const u128 min_include = static_cast<u128>(p) * p;
            if (static_cast<u128>(current_n) * min_include < static_cast<u128>(limit_exclusive)) {
                break;
            }
            ++i;
        }

        if (i >= npr) {
            add_solution(current_n);
            return;
        }

        const int p = primes_desc[static_cast<std::size_t>(i)];
        const int c = residue_mod3[static_cast<std::size_t>(p)];

        if (c == 0) {
            // Exclude p.
            dfs(i + 1, current_n);

            // Include p with exponent congruent to 2 (mod 3), minimal exponent 2.
            const int e0 = 2;
            const u64 power = pow_u64(static_cast<u64>(p), e0);
            if (static_cast<u128>(current_n) * static_cast<u128>(power) < static_cast<u128>(limit_exclusive)) {
                apply_factors(p, +1);
                included_primes.push_back(p);
                dfs(i + 1, current_n * power);
                included_primes.pop_back();
                apply_factors(p, -1);
            }
            return;
        }

        // Forced inclusion when residue is nonzero.
        const int m = (2 + c) % 3;
        const int e0 = (m == 0) ? 3 : m;
        const u64 power = pow_u64(static_cast<u64>(p), e0);
        if (static_cast<u128>(current_n) * static_cast<u128>(power) >= static_cast<u128>(limit_exclusive)) {
            return;
        }

        apply_factors(p, +1);
        included_primes.push_back(p);
        dfs(i + 1, current_n * power);
        included_primes.pop_back();
        apply_factors(p, -1);
    }

    u128 solve() {
        build_primes_and_spf();
        answer = 0;
        included_primes.clear();
        dfs(0, 1ULL);
        return answer;
    }
};

bool run_checkpoints() {
    // Problem statement example: n=50 qualifies.
    const u128 phi_50_sq = static_cast<u128>(50ULL) * 20ULL;  // phi(50)=20, and phi(n^2)=n*phi(n)
    if (!is_perfect_cube(phi_50_sq)) {
        std::cerr << "Checkpoint failed: n=50 should satisfy cube condition\n";
        return false;
    }

    // Cross-check on a smaller limit with brute force.
    const int small_limit = 200000;
    Solver solver_small(static_cast<u64>(small_limit));
    const u128 fast = solver_small.solve();
    const u128 slow = brute_sum(small_limit);
    if (fast != slow) {
        std::cerr << "Checkpoint failed: fast/brute mismatch for limit=" << small_limit << '\n';
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

    Solver solver(10000000000ULL);
    const u128 answer = solver.solve();
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
