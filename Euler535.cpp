#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'000ULL;

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) != 0 && (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

u128 sum_floor_sqrt(u64 n) {
    // Sum_{k=1..n} floor(sqrt(k)).
    // Count by threshold: floor(sqrt(k)) >= t iff k >= t^2.
    // => sum = sum_{t=1..s} (n - t^2 + 1) = s*(n+1) - sum_{t=1..s} t^2.
    const u64 s = isqrt_u64(n);
    const u128 su = static_cast<u128>(s);
    const u128 nn1 = static_cast<u128>(n) + 1;
    const u128 sumsq = su * (su + 1) * (2 * su + 1) / 6;
    return su * nn1 - sumsq;
}

u64 tri_mod(u64 n) {
    // n*(n+1)/2 mod 1e9, computed exactly in 128-bit.
    const u128 t = static_cast<u128>(n) * (static_cast<u128>(n) + 1) / 2;
    return static_cast<u64>(t % static_cast<u128>(kMod));
}

struct Solver {
    std::unordered_map<u64, u64> memoA;
    std::unordered_map<u64, u128> memoU;
    std::unordered_map<u64, u64> memoTmod;

    Solver() {
        memoA.reserve(1 << 20);
        memoU.reserve(1 << 20);
        memoTmod.reserve(1 << 20);
        memoA[0] = 0;
        memoU[0] = 0;
        memoTmod[0] = 0;
    }

    u64 A(u64 n) {
        auto it = memoA.find(n);
        if (it != memoA.end()) {
            return it->second;
        }
        if (n <= 1) {
            memoA.emplace(n, 0);
            return 0;
        }

        // A(n) is the largest k such that k + U(k) <= n.
        auto f = [&](u64 k) -> u128 { return static_cast<u128>(k) + U(k); };

        const u64 hi_max = n - 1;  // A(n) < n since S starts with a circled 1
        u64 lo = 0;
        u64 hi = 1;
        while (hi < hi_max && f(hi) <= static_cast<u128>(n)) {
            lo = hi;
            const u64 next = hi << 1;
            hi = (next < hi_max) ? next : hi_max;
        }
        if (f(hi) <= static_cast<u128>(n)) {
            memoA.emplace(n, hi);
            return hi;
        }
        while (lo + 1 < hi) {
            const u64 mid = lo + (hi - lo) / 2;
            if (f(mid) <= static_cast<u128>(n)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        memoA.emplace(n, lo);
        return lo;
    }

    u128 U(u64 n) {
        auto it = memoU.find(n);
        if (it != memoU.end()) {
            return it->second;
        }
        const u64 a = A(n);
        const u64 c = n - a;
        const u128 res = U(a) + sum_floor_sqrt(c);
        memoU.emplace(n, res);
        return res;
    }

    u64 T_mod(u64 n) {
        auto it = memoTmod.find(n);
        if (it != memoTmod.end()) {
            return it->second;
        }
        const u64 a = A(n);
        const u64 c = n - a;
        const u64 res = (T_mod(a) + tri_mod(c)) % kMod;
        memoTmod.emplace(n, res);
        return res;
    }

    u64 T_exact_u64(u64 n) {
        // Exact T(n) for moderate n where it fits in u64; used only for checkpoints.
        if (n == 0) {
            return 0;
        }
        const u64 a = A(n);
        const u64 c = n - a;
        const u128 res = static_cast<u128>(T_exact_u64(a)) +
                         static_cast<u128>(c) * (static_cast<u128>(c) + 1) / 2;
        return static_cast<u64>(res);
    }
};

bool run_checkpoints() {
    Solver s;
    if (s.T_exact_u64(1) != 1ULL) {
        std::cerr << "Checkpoint failed: T(1)\n";
        return false;
    }
    if (s.T_exact_u64(20) != 86ULL) {
        std::cerr << "Checkpoint failed: T(20) got " << s.T_exact_u64(20) << " A(20)=" << s.A(20)
                  << " T(9)=" << s.T_exact_u64(9) << " A(9)=" << s.A(9) << " T(4)=" << s.T_exact_u64(4)
                  << " A(4)=" << s.A(4) << " T(2)=" << s.T_exact_u64(2) << " A(2)=" << s.A(2) << '\n';
        return false;
    }
    if (s.T_exact_u64(1'000) != 364'089ULL) {
        std::cerr << "Checkpoint failed: T(10^3)\n";
        return false;
    }
    if (s.T_exact_u64(1'000'000'000ULL) != 498'676'527'978'348'241ULL) {
        std::cerr << "Checkpoint failed: T(10^9)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }
    Solver s;
    constexpr u64 n = 1'000'000'000'000'000'000ULL;
    const u64 ans = s.T_mod(n);
    std::cout << std::setw(9) << std::setfill('0') << ans << '\n';
    return 0;
}
