#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1'000'000'007LL;

int v2(const int x) {
    return __builtin_ctz(static_cast<unsigned>(x));
}

i64 mod_pow2(i64 exp) {
    i64 base = 2;
    i64 result = 1;
    while (exp > 0) {
        if (exp & 1LL) {
            result = static_cast<i64>((__int128)result * base % kMod);
        }
        base = static_cast<i64>((__int128)base * base % kMod);
        exp >>= 1LL;
    }
    return result;
}

i64 F(const int n, const int k) {
    const int g = std::gcd(n, k);
    if (v2(k) <= v2(n)) {
        return mod_pow2(static_cast<i64>(n - g + 1));
    }
    return mod_pow2(static_cast<i64>(n - g));
}

i64 solve(const int N) {
    std::vector<int> phi(static_cast<std::size_t>(N + 1), 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(N / 10));

    phi[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (phi[static_cast<std::size_t>(i)] == 0) {
            phi[static_cast<std::size_t>(i)] = i - 1;
            primes.push_back(i);
        }
        for (const int p : primes) {
            const i64 v = static_cast<i64>(i) * p;
            if (v > N) {
                break;
            }
            if (i % p == 0) {
                phi[static_cast<std::size_t>(v)] = phi[static_cast<std::size_t>(i)] * p;
                break;
            }
            phi[static_cast<std::size_t>(v)] = phi[static_cast<std::size_t>(i)] * (p - 1);
        }
    }

    std::vector<i64> pow2(static_cast<std::size_t>(N + 1), 1);
    for (int i = 1; i <= N; ++i) {
        pow2[static_cast<std::size_t>(i)] = (pow2[static_cast<std::size_t>(i - 1)] * 2) % kMod;
    }

    i64 ans = (2LL * N) % kMod;  // h = 1 contribution for each n

    for (int h = 2; h <= N; ++h) {
        i64 coeff;
        if ((h & 1) == 0) {
            coeff = (2LL * phi[static_cast<std::size_t>(h)]) % kMod;
        } else {
            coeff = (3LL * (phi[static_cast<std::size_t>(h)] / 2LL)) % kMod;
        }

        const int step = h - 1;
        const int count = N / h;
        int exp = step;
        for (int i = 0; i < count; ++i) {
            const i64 term = static_cast<i64>((__int128)coeff * pow2[static_cast<std::size_t>(exp)] % kMod);
            ans += term;
            if (ans >= kMod) {
                ans -= kMod;
            }
            exp += step;
        }
    }

    return ans;
}

}  // namespace

int main() {
    assert(F(3, 2) == 4);
    assert(F(8, 3) == 256);
    assert(F(9, 3) == 128);

    assert(solve(3) == 22);
    assert(solve(10) == 10'444);
    assert(solve(1'000) == 853'837'042);

    std::cout << solve(10'000'000) << '\n';
    return 0;
}
