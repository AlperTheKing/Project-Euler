#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

using u64 = std::uint64_t;

static long double harmonic_asymptotic(u64 n) {
    // Euler–Maclaurin: H_n = log(n) + gamma + 1/(2n) - 1/(12n^2) + O(1/n^4).
    static constexpr long double gamma = 0.577215664901532860606512090082402431L;
    const long double nn = static_cast<long double>(n);
    const long double inv = 1.0L / nn;
    const long double inv2 = inv * inv;
    return logl(nn) + gamma + 0.5L * inv - (1.0L / 12.0L) * inv2;
}

static long double B_large(u64 n) {
    // B_n = sum_{j=0..n-1} 2^{-j} / (n-j). The tail decays exponentially, so a small
    // fixed truncation gives far more than 8-decimal accuracy for n ~ 1e8.
    const int K = 200;
    long double sum = 0.0L;
    long double pow2 = 1.0L;  // 2^{-j}
    for (int j = 0; j <= K; ++j) {
        sum += pow2 / static_cast<long double>(n - static_cast<u64>(j));
        pow2 *= 0.5L;
    }
    return sum;
}

static long double S_small(int m) {
    // Exact O(m) for validation only.
    long double B = 0.0L;
    long double H = 0.0L;
    long double sum = 0.0L;
    for (int n = 1; n <= m; ++n) {
        H += 1.0L / static_cast<long double>(n);
        B = 1.0L / static_cast<long double>(n) + 0.5L * B;  // J_B(n)
        const long double JA = B - H / ldexpl(1.0L, n);      // H_n / 2^n
        const long double JB = B;
        sum += JA + JB;
    }
    return sum;
}

static long double S_large(u64 m) {
    // J_B(n) = (1/2^n) * sum_{k=1..n} 2^k/k, and satisfies B_n = 1/n + B_{n-1}/2.
    // J_A(n) = J_B(n) - H_n/2^n, hence:
    //   J_A(n)+J_B(n) = 2B_n - H_n/2^n.
    //
    // Summing:
    //   sum_{n<=m} B_n = 2 H_m - B_m
    // and sum_{n>=1} H_n/2^n = 2 ln 2, with the tail beyond m being ~ H_m/2^m (negligible here).
    const long double Hm = harmonic_asymptotic(m);
    const long double Bm = B_large(m);
    const long double two_ln2 = 2.0L * logl(2.0L);
    return 4.0L * Hm - 2.0L * Bm - two_ln2;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Statement validation at n=6 (8 decimals).
    {
        long double B = 0.0L;
        long double H = 0.0L;
        for (int n = 1; n <= 6; ++n) {
            H += 1.0L / static_cast<long double>(n);
            B = 1.0L / static_cast<long double>(n) + 0.5L * B;
        }
        const long double JB6 = B;
        const long double JA6 = B - H / ldexpl(1.0L, 6);
        auto r8 = [](long double x) { return std::llround(x * 100000000.0L); };
        assert(r8(JA6) == 39505208);
        assert(r8(JB6) == 43333333);
        assert(r8(S_small(6)) == 758932292);
    }
    // Sanity: for moderately large m, the tail of sum H_n/2^n is negligible, so S_large(m) matches the exact DP.
    {
        const int m = 1000;
        const long double exact = S_small(m);
        const long double approx = S_large(static_cast<u64>(m));
        assert(fabsl(exact - approx) < 1e-12L);
    }

    const u64 m = 123456789ULL;
    const long double ans = S_large(m);
    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(ans) << '\n';
    return 0;
}
