#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

constexpr std::uint64_t kMod = 1'000'000'007ULL;

inline int parity32(std::uint32_t x) {
    return __builtin_parity(x);
}

std::vector<std::uint32_t> build_spf(int n) {
    std::vector<std::uint32_t> spf(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::uint32_t> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));

    spf[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = static_cast<std::uint32_t>(i);
            primes.push_back(static_cast<std::uint32_t>(i));
        }
        for (std::uint32_t p : primes) {
            const std::uint64_t v = static_cast<std::uint64_t>(i) * p;
            if (v > static_cast<std::uint64_t>(n)) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
            if (p == spf[static_cast<std::size_t>(i)]) {
                break;
            }
        }
    }
    return spf;
}

std::uint64_t solve(int n) {
    const std::vector<std::uint32_t> spf = build_spf(n);
    std::vector<std::uint32_t> inert_exp(static_cast<std::size_t>(n / 2) + 1, 0);

    std::uint32_t exp2 = 0;
    int parity = 0;

    std::uint64_t fact_mod = 1;
    std::uint64_t sum_mod = 0;

    for (int k = 1; k <= n; ++k) {
        fact_mod = (fact_mod * static_cast<std::uint64_t>(k)) % kMod;

        int x = k;
        while (x > 1) {
            const std::uint32_t p = spf[static_cast<std::size_t>(x)];
            std::uint32_t cnt = 0;
            do {
                x /= static_cast<int>(p);
                ++cnt;
            } while (x > 1 && spf[static_cast<std::size_t>(x)] == p);

            if (p == 2U) {
                const int old_bit = parity32(exp2);
                exp2 += cnt;
                const int new_bit = parity32(exp2);
                parity ^= (old_bit ^ new_bit);
            } else {
                const int r = static_cast<int>(p & 7U);
                if (r == 5 || r == 7) {
                    const std::size_t idx = static_cast<std::size_t>(p >> 1U);
                    const std::uint32_t old_exp = inert_exp[idx];
                    const int old_bit = parity32(old_exp);
                    const std::uint32_t new_exp = old_exp + cnt;
                    inert_exp[idx] = new_exp;
                    const int new_bit = parity32(new_exp);
                    parity ^= (old_bit ^ new_bit);
                }
            }
        }

        if (parity == 0) {
            sum_mod += fact_mod;
            if (sum_mod >= kMod) {
                sum_mod -= kMod;
            }
        }
    }

    return sum_mod;
}

void run_validations() {
    assert(solve(4) == 25ULL);
    assert(solve(7) == 745ULL);
    assert(solve(100) == 709'772'949ULL);
}

}  // namespace

int main() {
    run_validations();
    constexpr int kN = 100'000'000;
    std::cout << solve(kN) << '\n';
    return 0;
}
