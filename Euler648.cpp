#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr int DEG = 1000;
constexpr u64 kMod = 1'000'000'000ULL;

inline std::size_t idxC(const int n, const int k) { return static_cast<std::size_t>(n) * (DEG + 1) + k; }

std::vector<u32> binom_table() {
    const int MAXN = 2 * DEG;
    std::vector<u32> C(static_cast<std::size_t>(MAXN + 1) * (DEG + 1), 0U);
    C[idxC(0, 0)] = 1U;
    for (int n = 1; n <= MAXN; ++n) {
        C[idxC(n, 0)] = 1U;
        const int up = std::min(n, DEG);
        for (int k = 1; k <= up; ++k) {
            const u64 a = C[idxC(n - 1, k)];
            const u64 b = C[idxC(n - 1, k - 1)];
            C[idxC(n, k)] = static_cast<u32>((a + b) % kMod);
        }
    }
    return C;
}

u64 solve() {
    const auto C = binom_table();

    std::vector<u64> S(DEG + 1, 0), g(DEG + 1, 0), P(DEG + 1, 0), y(DEG + 1, 0), newP(DEG + 1, 0);
    std::vector<u128> acc(DEG + 1);

    P[0] = 1;
    y[0] = 1;

    for (int n = 1; n <= DEG; ++n) {
        const int m = 2 * (n - 1);
        const int up = std::min(m, DEG);
        for (int k = 0; k <= up; ++k) {
            u64 coef = C[idxC(m, k)];
            if (k & 1) {
                coef = (coef == 0) ? 0 : (kMod - coef);
            }
            S[k] += coef;
            if (S[k] >= kMod) S[k] -= kMod;
        }

        std::fill(g.begin(), g.end(), 0);
        for (int k = 1; k <= DEG; ++k) g[k] = S[k - 1];
        for (int k = 2; k <= DEG; ++k) {
            g[k] += kMod - S[k - 2];
            if (g[k] >= kMod) g[k] -= kMod;
        }

        std::fill(acc.begin(), acc.end(), 0);
        for (int i = 0; i <= DEG; ++i) {
            const u64 Pi = P[i];
            if (Pi == 0) continue;
            const int maxj = DEG - i;
            for (int j = 1; j <= maxj; ++j) {
                const u64 gj = g[j];
                if (gj == 0) continue;
                acc[i + j] += (u128)Pi * (u128)gj;
            }
        }

        for (int k = 0; k <= DEG; ++k) newP[k] = static_cast<u64>(acc[k] % (u128)kMod);
        P.swap(newP);

        for (int k = 0; k <= DEG; ++k) {
            y[k] += P[k];
            y[k] %= kMod;
        }
    }

    assert(y[0] == 1);
    assert(y[1] == 1);
    assert(y[10] == 53964);
    assert(y[50] == 842'418'857ULL);
    assert((y[5] + kMod - y[4]) % kMod == kMod - 18);
    assert((y[10] + kMod - y[9]) % kMod == 45'176ULL);

    return y[DEG] % kMod;
}

}  // namespace

int main() {
    std::cout << solve() << "\n";
    return 0;
}
