#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using i128 = __int128_t;

std::int64_t lcm_i64(const std::int64_t a, const std::int64_t b) {
    return a / std::gcd(a, b) * b;
}

// Compute E(n) exactly via a common denominator D = lcm(1..N).
// The derived recurrence is:
//   E(0)=0,
//   E(n) = (sum_{k=0..n-1} E(k) + (2^n - n - 1)) / n.
struct ExpectedMoves {
    int N = 0;
    std::int64_t D = 1;          // common denominator
    std::vector<i128> scaled_E;  // E(n) * D

    explicit ExpectedMoves(const int N_) : N(N_) {
        for (int i = 1; i <= N; ++i) {
            D = lcm_i64(D, i);
        }
        scaled_E.assign(static_cast<std::size_t>(N + 1), 0);

        i128 prefix = 0;  // sum_{k=0..n-1} E(k) * D
        for (int n = 1; n <= N; ++n) {
            const std::int64_t two_pow_n = (n < 63) ? (1LL << n) : 0LL;
            assert(n < 63);  // sufficient for N=30 here
            const std::int64_t c = two_pow_n - n - 1;  // (2^n - n - 1)

            const i128 numerator = prefix + static_cast<i128>(c) * static_cast<i128>(D);
            assert(numerator % n == 0);
            scaled_E[static_cast<std::size_t>(n)] = numerator / n;
            prefix += scaled_E[static_cast<std::size_t>(n)];
        }
    }

    i128 scaled(int n) const { return scaled_E[static_cast<std::size_t>(n)]; }
};

bool run_checkpoints() {
    const ExpectedMoves em(30);
    const i128 D = static_cast<i128>(em.D);

    // Given checkpoints from the statement.
    {
        const i128 expect = (D / 4) * 13;  // 13/4
        if (em.scaled(4) != expect) {
            std::cerr << "Checkpoint failed: E(4)\n";
            return false;
        }
    }
    {
        const i128 expect = (D / 40) * 4629;  // 4629/40
        if (em.scaled(10) != expect) {
            std::cerr << "Checkpoint failed: E(10)\n";
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    const ExpectedMoves em(30);
    const i128 D = static_cast<i128>(em.D);
    const i128 num100 = em.scaled(30) * 100;

    i128 q = num100 / D;
    const i128 rem = num100 % D;
    if (rem * 2 >= D) {
        ++q;
    }

    const std::int64_t q64 = static_cast<std::int64_t>(q);
    const std::int64_t integer = q64 / 100;
    const std::int64_t frac = q64 % 100;
    std::cout << integer << '.' << std::setw(2) << std::setfill('0') << frac << '\n';
    return 0;
}

