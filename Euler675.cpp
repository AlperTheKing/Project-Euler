#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr i64 MOD = 1'000'000'087LL;

std::vector<int> build_spf(int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));

    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(p) * static_cast<i64>(i);
            if (v > n || p > spf[static_cast<std::size_t>(i)]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }

    return spf;
}

int max_term_bound(int n) {
    int e2 = 0;
    for (int x = n; x > 0; x /= 2) {
        e2 += x / 2;
    }
    return 2 * e2 + 1;
}

std::vector<i64> build_inverses(int bound) {
    std::vector<i64> inv(static_cast<std::size_t>(bound + 1), 0LL);
    inv[1] = 1LL;
    for (int i = 2; i <= bound; ++i) {
        inv[static_cast<std::size_t>(i)] =
            MOD - static_cast<i64>((static_cast<__int128>(MOD / i) *
                                    inv[static_cast<std::size_t>(MOD % i)]) %
                                   MOD);
    }
    return inv;
}

i64 S_of_number(u64 n) {
    if (n == 1ULL) {
        return 1LL;
    }
    i64 out = 1LL;
    u64 x = n;
    for (u64 p = 2ULL; p * p <= x; ++p) {
        if (x % p != 0ULL) {
            continue;
        }
        int e = 0;
        while (x % p == 0ULL) {
            x /= p;
            ++e;
        }
        out = static_cast<i64>((static_cast<__int128>(out) * (1 + 2 * e)) % MOD);
    }
    if (x > 1ULL) {
        out = static_cast<i64>((static_cast<__int128>(out) * 3LL) % MOD);
    }
    return out;
}

i64 solve(int n) {
    const std::vector<int> spf = build_spf(n);
    const int bound = max_term_bound(n);
    const std::vector<i64> inv = build_inverses(bound);

    std::vector<int> exp(static_cast<std::size_t>(n + 1), 0);

    i64 current = 1LL;
    i64 ans = 0LL;

    for (int i = 2; i <= n; ++i) {
        int x = i;
        while (x > 1) {
            const int p = spf[static_cast<std::size_t>(x)];
            int c = 0;
            do {
                x /= p;
                ++c;
            } while (x % p == 0);

            const int old_term = 1 + 2 * exp[static_cast<std::size_t>(p)];
            exp[static_cast<std::size_t>(p)] += c;
            const int new_term = 1 + 2 * exp[static_cast<std::size_t>(p)];

            current = static_cast<i64>((static_cast<__int128>(current) *
                                        static_cast<i64>(new_term)) %
                                       MOD);
            current = static_cast<i64>((static_cast<__int128>(current) *
                                        inv[static_cast<std::size_t>(old_term)]) %
                                       MOD);
        }

        ans += current;
        if (ans >= MOD) {
            ans -= MOD;
        }
    }

    return ans;
}

}  // namespace

int main() {
    assert(S_of_number(6ULL) == 9LL);
    assert(solve(10) == 4'821LL);

    std::cout << solve(10'000'000) << "\n";
    return 0;
}
