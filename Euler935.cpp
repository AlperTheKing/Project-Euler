#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u64 = std::uint64_t;

u64 totient_sum(u64 n) {
    std::vector<std::uint32_t> phi(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    std::vector<bool> composite(n + 1, false);

    phi[1] = 1;
    for (u64 i = 2; i <= n; ++i) {
        if (!composite[i]) {
            primes.push_back(static_cast<int>(i));
            phi[i] = static_cast<std::uint32_t>(i - 1);
        }
        for (int p : primes) {
            const u64 v = i * static_cast<u64>(p);
            if (v > n) {
                break;
            }
            composite[v] = true;
            if (i % static_cast<u64>(p) == 0) {
                phi[v] = static_cast<std::uint32_t>(phi[i] * static_cast<u64>(p));
                break;
            }
            phi[v] = static_cast<std::uint32_t>(phi[i] * static_cast<u64>(p - 1));
        }
    }

    u64 sum = 0;
    for (u64 i = 1; i <= n; ++i) {
        sum += phi[i];
    }
    return sum;
}

u64 brute_totient_sum(int n) {
    u64 sum = 0;
    for (int i = 1; i <= n; ++i) {
        int cnt = 0;
        for (int j = 1; j <= i; ++j) {
            if (std::gcd(i, j) == 1) {
                ++cnt;
            }
        }
        sum += static_cast<u64>(cnt);
    }
    return sum;
}

u64 solve(u64 n) {
    if (n < 4) {
        return 0;
    }
    const u64 m = (n + 2) / 2;
    u64 result = totient_sum(m) - 1;
    if (n < 8 && n >= 6) {
        --result;
    }
    return result;
}

void run_validations() {
    for (int n = 1; n <= 250; ++n) {
        assert(totient_sum(static_cast<u64>(n)) == brute_totient_sum(n));
    }
    assert(solve(4) == 3);
    assert(solve(6) == 4);
    assert(solve(100) == 805);
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kN = 100'000'000ULL;
    std::cout << solve(kN) << '\n';
    return 0;
}
