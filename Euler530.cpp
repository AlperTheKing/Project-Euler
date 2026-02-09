#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using u32 = std::uint32_t;

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) != 0 && (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

u64 icbrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while ((r + 1) != 0 && (r + 1) * (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r * r > x) {
        --r;
    }
    return r;
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const u64 digit = static_cast<u64>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

// D(n) = sum_{m<=n} tau(m) = sum_{i=1..n} floor(n/i).
u64 divisor_summatory(const u64 n) {
    u128 res = 0;
    for (u64 l = 1; l <= n;) {
        const u64 q = n / l;
        const u64 r = n / q;
        res += static_cast<u128>(q) * static_cast<u128>(r - l + 1);
        l = r + 1;
    }
    return static_cast<u64>(res);
}

std::vector<u32> totient_sieve(const u64 n) {
    std::vector<u32> phi(static_cast<std::size_t>(n + 1), 0);
    for (u64 i = 0; i <= n; ++i) {
        phi[static_cast<std::size_t>(i)] = static_cast<u32>(i);
    }
    if (n >= 1) {
        phi[1] = 1;
    }
    for (u64 p = 2; p <= n; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (u64 m = p; m <= n; m += p) {
            const u32 cur = phi[static_cast<std::size_t>(m)];
            phi[static_cast<std::size_t>(m)] = static_cast<u32>(cur - cur / static_cast<u32>(p));
        }
    }
    return phi;
}

u128 compute_F(const u64 N) {
    const u64 limit = isqrt_u64(N);
    const u64 K = icbrt_u64(N);

    const auto phi = totient_sieve(limit);

    const u64 max_small = N / (K * K);  // ~= K
    std::vector<u32> tau(static_cast<std::size_t>(max_small + 1), 0);
    for (u64 d = 1; d <= max_small; ++d) {
        for (u64 m = d; m <= max_small; m += d) {
            ++tau[static_cast<std::size_t>(m)];
        }
    }
    std::vector<u64> D_small(static_cast<std::size_t>(max_small + 1), 0);
    for (u64 i = 1; i <= max_small; ++i) {
        D_small[static_cast<std::size_t>(i)] =
            D_small[static_cast<std::size_t>(i - 1)] + static_cast<u64>(tau[static_cast<std::size_t>(i)]);
    }

    // F(N) = sum_{t<=sqrt(N)} phi(t) * D(floor(N/t^2)).
    u128 ans = 0;

    // For t <= K, floor(N/t^2) is large; compute D via divisor_summatory.
    // Group equal quotients to avoid recomputing D for duplicates.
    u64 t = 1;
    while (t <= K) {
        const u64 v = N / (t * t);
        u64 r = isqrt_u64(N / v);
        if (r > K) {
            r = K;
        }
        u64 sum_phi = 0;
        for (u64 i = t; i <= r; ++i) {
            sum_phi += static_cast<u64>(phi[static_cast<std::size_t>(i)]);
        }
        const u64 Dv = divisor_summatory(v);
        ans += static_cast<u128>(sum_phi) * static_cast<u128>(Dv);
        t = r + 1;
    }

    // For t > K, v <= max_small, so D can be tabled.
    for (u64 i = K + 1; i <= limit; ++i) {
        const u64 v = N / (i * i);
        ans += static_cast<u128>(phi[static_cast<std::size_t>(i)]) *
               static_cast<u128>(D_small[static_cast<std::size_t>(v)]);
    }

    return ans;
}

bool run_checkpoints() {
    if (compute_F(10) != static_cast<u128>(32)) {
        std::cerr << "Checkpoint failed: F(10)\n";
        return false;
    }
    if (compute_F(1000) != static_cast<u128>(12776)) {
        std::cerr << "Checkpoint failed: F(1000)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 N = 1'000'000'000'000'000ULL;
    const u128 ans = compute_F(N);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}

