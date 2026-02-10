#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kInv2 = 500'000'004ULL;

struct PhiSummatory {
    int limit = 0;
    std::vector<int> primes;
    std::vector<int> lp;
    std::vector<u32> phi;
    std::vector<u64> pref;
    std::unordered_map<u64, u64> memo;

    explicit PhiSummatory(const int n) : limit(n), lp(n + 1, 0), phi(n + 1, 0), pref(n + 1, 0) {
        primes.reserve(static_cast<std::size_t>(n / 10));
        phi[1] = 1;
        for (int i = 2; i <= n; ++i) {
            if (lp[i] == 0) {
                lp[i] = i;
                primes.push_back(i);
                phi[i] = static_cast<u32>(i - 1);
            }
            for (int p : primes) {
                if (p > lp[i] || (u64)i * (u64)p > (u64)n) break;
                lp[i * p] = p;
                if (p == lp[i]) {
                    phi[i * p] = phi[i] * static_cast<u32>(p);
                    break;
                }
                phi[i * p] = phi[i] * static_cast<u32>(p - 1);
            }
        }
        for (int i = 1; i <= n; ++i) {
            pref[i] = (pref[i - 1] + static_cast<u64>(phi[i])) % kMod;
        }
        memo.reserve(1 << 20);
    }

    u64 sum_phi(const u64 n) {
        if (n <= static_cast<u64>(limit)) return pref[static_cast<std::size_t>(n)];
        const auto it = memo.find(n);
        if (it != memo.end()) return it->second;

        const u64 nn = n % kMod;
        u64 res = nn * ((n + 1) % kMod) % kMod;
        res = res * kInv2 % kMod;

        for (u64 l = 2; l <= n;) {
            const u64 q = n / l;
            const u64 r = n / q;
            const u64 cnt = (r - l + 1) % kMod;
            res = (res + kMod - cnt * sum_phi(q) % kMod) % kMod;
            l = r + 1;
        }

        memo.emplace(n, res);
        return res;
    }
};

u64 solve(const u64 n, PhiSummatory& ph) {
    u64 ans = 0;
    for (u64 m = n / 2; m >= 2; m /= 2) {
        const u64 add = (ph.sum_phi(m) + kMod - 1) % kMod;
        ans += add;
        ans %= kMod;
    }
    return ans;
}

u64 brute(const u32 n) {
    u64 cnt = 0;
    for (u32 a = 1; a <= n; ++a) {
        for (u32 b = a + 1; b <= n; ++b) {
            const u32 g = std::gcd(a, b);
            if (g > 1 && (g & (g - 1)) == 0) ++cnt;
        }
    }
    return cnt % kMod;
}

}  // namespace

int main() {
    PhiSummatory ph(5'000'000);

    assert(solve(100, ph) == 1031);
    assert(solve(1'000'000, ph) == 321'418'433ULL);
    assert(solve(2000, ph) == brute(2000));

    std::cout << solve(100'000'000'000ULL, ph) << "\n";
    return 0;
}
