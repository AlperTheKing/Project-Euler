#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

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
            const u64 v = static_cast<u64>(p) * static_cast<u64>(i);
            if (v > static_cast<u64>(n) || p > spf[static_cast<std::size_t>(i)]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }

    return spf;
}

std::vector<std::pair<u64, int>> factor_square(int m, const std::vector<int>& spf) {
    std::vector<std::pair<u64, int>> fac;
    int x = m;
    while (x > 1) {
        const int p = spf[static_cast<std::size_t>(x)];
        int e = 0;
        do {
            x /= p;
            ++e;
        } while (x % p == 0);
        fac.push_back({static_cast<u64>(p), 2 * e});
    }
    return fac;
}

void build_divisors_rec(int idx, u64 cur, const std::vector<std::pair<u64, int>>& fac,
                        std::vector<u64>& out) {
    if (idx == static_cast<int>(fac.size())) {
        out.push_back(cur);
        return;
    }

    const u64 p = fac[static_cast<std::size_t>(idx)].first;
    const int e = fac[static_cast<std::size_t>(idx)].second;
    u64 mul = 1ULL;
    for (int i = 0; i <= e; ++i) {
        build_divisors_rec(idx + 1, cur * mul, fac, out);
        mul *= p;
    }
}

u64 SP(int n) {
    const std::vector<int> spf = build_spf(n);

    u128 ans = 0;
    std::vector<u64> divs;

    for (int m = 1; m <= n; ++m) {
        const u64 n2 = static_cast<u64>(m) * static_cast<u64>(m);

        const auto fac = factor_square(m, spf);
        divs.clear();
        divs.reserve(256);
        build_divisors_rec(0, 1ULL, fac, divs);
        std::sort(divs.begin(), divs.end());

        const int L = static_cast<int>(divs.size());

        for (int i = 0; i < L; ++i) {
            const u64 w = divs[static_cast<std::size_t>(i)];

            for (int j = i; j < L; ++j) {
                const u64 z = divs[static_cast<std::size_t>(j)];
                const u64 wz = w * z;
                if (wz > n2) {
                    break;
                }
                if (n2 % wz != 0ULL) {
                    continue;
                }

                const u64 rem = n2 / wz;
                if (z * z > rem) {
                    break;
                }

                for (int k = j; k < L; ++k) {
                    const u64 y = divs[static_cast<std::size_t>(k)];
                    if (y * y > rem) {
                        break;
                    }
                    if (rem % y != 0ULL) {
                        continue;
                    }

                    const u64 x = rem / y;
                    if (x < y) {
                        continue;
                    }
                    if (x >= y + z + w) {
                        continue;
                    }
                    const u64 p = x + y + z + w;
                    if ((p & 1ULL) != 0ULL) {
                        continue;
                    }

                    ans += static_cast<u128>(p);
                }
            }
        }
    }

    return static_cast<u64>(ans);
}

}  // namespace

int main() {
    assert(SP(10) == 186ULL);
    assert(SP(100) == 23'238ULL);

    std::cout << SP(1'000'000) << "\n";
    return 0;
}
