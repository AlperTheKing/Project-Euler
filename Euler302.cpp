#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;

constexpr int kPrimeSieveN = 6 * 105000;
int prime_count = 0;
int primes[kPrimeSieveN / 4];
unsigned char erat[kPrimeSieveN + 1];

int divisors[kPrimeSieveN + 1][24];
unsigned char divisors_count[kPrimeSieveN + 1];

u64 answer = 0;

int gcd_small(int a, int b) {
    while (b != 0) {
        const int c = a % b;
        a = b;
        b = c;
    }
    return a;
}

int cube_root_floor(const u64 n) {
    if (n <= 1ULL) {
        return static_cast<int>(n);
    }
    int r = static_cast<int>(std::pow(static_cast<long double>(n), 1.0L / 3.0L));
    while (static_cast<u64>(r + 1) * static_cast<u64>(r + 1) * static_cast<u64>(r + 1) <= n) ++r;
    while (static_cast<u64>(r) * static_cast<u64>(r) * static_cast<u64>(r) > n) --r;
    return r;
}

void count_recursive(const int max_prime_idx,
                     const int* active_factors,
                     const int active_count,
                     int gcd_exp_n,
                     int gcd_exp_phi,
                     const u64 lim) {
    if (lim == 0ULL) {
        return;
    }

    if (gcd_exp_n == 1) {
        int g = gcd_exp_phi;
        int i = 0;
        while (i < active_count) {
            int j = i + 1;
            while (j < active_count && active_factors[j] == active_factors[i]) ++j;
            const int cnt = j - i;
            if (cnt < 2) {
                g = 0;
                break;
            }
            g = gcd_small(g, cnt);
            i = j;
        }
        if (g == 1) {
            ++answer;
        }
    }

    int min_idx = 0;
    if (active_count >= 2) {
        if (active_factors[active_count - 1] != active_factors[active_count - 2]) {
            min_idx = active_factors[active_count - 1];
        } else {
            for (int k = active_count - 2; k >= 1; --k) {
                if (active_factors[k] != active_factors[k + 1] && active_factors[k] != active_factors[k - 1]) {
                    min_idx = active_factors[k];
                    break;
                }
            }
        }
    } else if (active_count == 1) {
        min_idx = active_factors[0];
    }

    int max_idx = cube_root_floor(lim);
    if (primes[max_prime_idx - 1] > max_idx) {
        int a = 0;
        int b = prime_count - 1;
        while (a < b - 1) {
            const int c = (a + b) / 2;
            if (primes[c] > max_idx) {
                b = c;
            } else {
                a = c;
            }
        }
        max_idx = a;
    }
    if (max_idx > max_prime_idx - 1) {
        max_idx = max_prime_idx - 1;
    }

    int merged[70];
    for (int j = active_count - 1; j >= 0;) {
        int overdeg = 0;
        for (int y = j; y >= 0 && active_factors[y] == active_factors[j]; --y) ++overdeg;

        if (active_factors[j] < max_prime_idx) {
            int k1 = 0;
            int k2 = 0;
            int merged_count = 0;
            const int idx = active_factors[j];
            const int p_minus_one = primes[idx] - 1;
            const int div_cnt = divisors_count[p_minus_one];

            while (k1 < active_count && k2 < div_cnt) {
                if (active_factors[k1] < divisors[p_minus_one][k2]) {
                    merged[merged_count++] = active_factors[k1++];
                } else {
                    merged[merged_count++] = divisors[p_minus_one][k2++];
                }
            }
            while (k2 < div_cnt) merged[merged_count++] = divisors[p_minus_one][k2++];
            while (k1 < active_count) merged[merged_count++] = active_factors[k1++];

            int g2 = gcd_exp_phi;
            while (merged_count > 0 && merged[merged_count - 1] > active_factors[j]) {
                int y = merged_count - 1;
                while (y >= 0 && merged[y] == merged[merged_count - 1]) --y;
                g2 = gcd_small(g2, merged_count - y - 1);
                merged_count = y + 1;
            }
            if (merged_count > 0 && merged[merged_count - 1] == active_factors[j]) {
                merged_count -= overdeg;
            }

            const u64 p = static_cast<u64>(primes[idx]);
            u64 pow = p;
            int deg = 2;
            while (pow <= lim / p) {
                const int new_g1 = gcd_small(gcd_exp_n, deg);
                const int new_g2 = gcd_small(g2, deg - 1 + overdeg);
                count_recursive(idx, merged, merged_count, new_g1, new_g2, lim / pow / p);

                if (pow > (std::numeric_limits<u64>::max() / p)) break;
                pow *= p;
                ++deg;
            }
        }

        j -= overdeg;
        if (overdeg == 1) break;
    }

    int now = 0;
    while (now < active_count && active_factors[now] < min_idx) ++now;

    for (int idx = min_idx; idx <= max_idx; ++idx) {
        if (now < active_count && idx == active_factors[now]) {
            while (now < active_count && active_factors[now] == idx) ++now;
            continue;
        }

        int k1 = 0;
        int k2 = 0;
        int merged_count = 0;
        const int p_minus_one = primes[idx] - 1;
        const int div_cnt = divisors_count[p_minus_one];
        while (k1 < active_count && k2 < div_cnt) {
            if (active_factors[k1] < divisors[p_minus_one][k2]) {
                merged[merged_count++] = active_factors[k1++];
            } else {
                merged[merged_count++] = divisors[p_minus_one][k2++];
            }
        }
        while (k2 < div_cnt) merged[merged_count++] = divisors[p_minus_one][k2++];
        while (k1 < active_count) merged[merged_count++] = active_factors[k1++];

        int g2 = gcd_exp_phi;
        while (merged_count > 0 && merged[merged_count - 1] > idx) {
            int y = merged_count - 1;
            while (y >= 0 && merged[y] == merged[merged_count - 1]) --y;
            g2 = gcd_small(g2, merged_count - y - 1);
            merged_count = y + 1;
        }

        const u64 p = static_cast<u64>(primes[idx]);
        u64 pow = p * p;
        int deg = 3;
        while (pow <= lim / p) {
            const int new_g1 = gcd_small(gcd_exp_n, deg);
            const int new_g2 = gcd_small(g2, deg - 1);
            count_recursive(idx, merged, merged_count, new_g1, new_g2, lim / pow / p);

            if (pow > (std::numeric_limits<u64>::max() / p)) break;
            pow *= p;
            ++deg;
        }
    }
}

void prepare_data() {
    primes[0] = 2;
    primes[1] = 3;
    prime_count = 2;
    erat[1] = 1;

    const int sq = static_cast<int>(std::sqrt(static_cast<long double>(kPrimeSieveN)));
    for (int i = 0; i <= sq; i += 6) {
        if (erat[i + 1] == 0) {
            const int p = i + 1;
            const int step = p * 6;
            for (int j = p * p; j <= kPrimeSieveN; j += step) erat[j] = 1;
            for (int j = p * (p + 4); j <= kPrimeSieveN; j += step) erat[j] = 1;
            primes[prime_count++] = p;
        }
        if (erat[i + 5] == 0) {
            const int p = i + 5;
            const int step = p * 6;
            for (int j = p * p; j <= kPrimeSieveN; j += step) erat[j] = 1;
            for (int j = p * (p + 2); j <= kPrimeSieveN; j += step) erat[j] = 1;
            primes[prime_count++] = p;
        }
    }

    const int start = sq - (sq % 6) + 6;
    for (int i = start; i < kPrimeSieveN; i += 6) {
        if (erat[i + 1] == 0) primes[prime_count++] = i + 1;
        if (erat[i + 5] == 0) primes[prime_count++] = i + 5;
    }

    for (int i = 0; i < prime_count; ++i) {
        for (u64 deg = static_cast<u64>(primes[i]); deg <= static_cast<u64>(kPrimeSieveN);) {
            for (int j = static_cast<int>(deg); j <= kPrimeSieveN; j += static_cast<int>(deg)) {
                divisors[j][divisors_count[j]++] = i;
            }
            if (deg <= static_cast<u64>(kPrimeSieveN / primes[i])) {
                deg *= static_cast<u64>(primes[i]);
            } else {
                break;
            }
        }
    }
}

u64 solve(const u64 limit) {
    answer = 0;
    int root_factors[70] = {0};
    count_recursive(prime_count, root_factors, 0, 0, 0, limit);
    return answer;
}

}  // namespace

int main() {
    prepare_data();
    assert(solve(10'000ULL) == 7ULL);
    assert(solve(100'000'000ULL) == 656ULL);
    std::cout << solve(1'000'000'000'000'000'000ULL) << '\n';
    return 0;
}
