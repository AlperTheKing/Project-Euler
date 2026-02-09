#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

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

// Extended gcd: returns (g,x,y) with ax+by=g.
i64 ext_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0, y1 = 0;
    const i64 g = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

// Smallest non-negative solution to:
//   x ≡ a (mod n)
//   x ≡ b (mod m)
// or 0 if no solution.
u64 crt_min_nonneg(u64 a, u64 n, u64 b, u64 m) {
    const u64 g = std::gcd(n, m);
    const i64 diff = static_cast<i64>(b) - static_cast<i64>(a);
    if (diff % static_cast<i64>(g) != 0) {
        return 0;
    }

    const u64 n1 = n / g;
    const u64 m1 = m / g;
    const i64 rhs = diff / static_cast<i64>(g);

    // Solve n1 * t ≡ rhs (mod m1), with gcd(n1,m1)=1.
    i64 x = 0, y = 0;
    const i64 gg = ext_gcd(static_cast<i64>(n1), static_cast<i64>(m1), x, y);
    (void)gg;

    i64 inv = x % static_cast<i64>(m1);
    if (inv < 0) {
        inv += static_cast<i64>(m1);
    }
    i64 rhs_mod = rhs % static_cast<i64>(m1);
    if (rhs_mod < 0) {
        rhs_mod += static_cast<i64>(m1);
    }
    const u64 t = static_cast<u64>(static_cast<u128>(rhs_mod) * static_cast<u64>(inv) % m1);

    // Minimal solution is x = a + n * t in [0, lcm).
    return a + n * t;
}

std::vector<u64> compute_phi_up_to(const int N) {
    std::vector<u64> phi(static_cast<std::size_t>(N + 1), 0);
    for (int i = 0; i <= N; ++i) {
        phi[static_cast<std::size_t>(i)] = static_cast<u64>(i);
    }
    if (N >= 1) {
        phi[1] = 1;
    }
    for (int p = 2; p <= N; ++p) {
        if (phi[static_cast<std::size_t>(p)] != static_cast<u64>(p)) {
            continue;
        }
        for (int m = p; m <= N; m += p) {
            phi[static_cast<std::size_t>(m)] -= phi[static_cast<std::size_t>(m)] / static_cast<u64>(p);
        }
    }
    return phi;
}

u128 solve() {
    constexpr int L = 1'000'000;
    constexpr int R_exclusive = 1'005'000;
    constexpr int N = R_exclusive - L;  // 5000

    const auto phi_all = compute_phi_up_to(R_exclusive);

    std::vector<u64> nvals(static_cast<std::size_t>(N));
    std::vector<u64> phis(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        const int n = L + i;
        nvals[static_cast<std::size_t>(i)] = static_cast<u64>(n);
        phis[static_cast<std::size_t>(i)] = phi_all[static_cast<std::size_t>(n)];
    }

    u128 sum = 0;
    for (int i = 0; i < N; ++i) {
        const u64 n = nvals[static_cast<std::size_t>(i)];
        const u64 a = phis[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < N; ++j) {
            const u64 m = nvals[static_cast<std::size_t>(j)];
            const u64 b = phis[static_cast<std::size_t>(j)];
            const u64 x = crt_min_nonneg(a, n, b, m);
            sum += static_cast<u128>(x);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (crt_min_nonneg(2, 4, 4, 6) != 10) {
        std::cerr << "Checkpoint failed: g(2,4,4,6)\n";
        return false;
    }
    if (crt_min_nonneg(3, 4, 4, 6) != 0) {
        std::cerr << "Checkpoint failed: g(3,4,4,6)\n";
        return false;
    }
    // f(4,6) = g(phi(4),4,phi(6),6) = g(2,4,2,6) = 2.
    if (crt_min_nonneg(2, 4, 2, 6) != 2) {
        std::cerr << "Checkpoint failed: f(4,6)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }
    std::cout << to_string_u128(solve()) << '\n';
    return 0;
}

