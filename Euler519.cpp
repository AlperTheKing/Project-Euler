#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'000'000'000ULL;  // last 9 digits

u64 mod_norm(i128 v) {
    const i128 m = static_cast<i128>(kMod);
    v %= m;
    if (v < 0) {
        v += m;
    }
    return static_cast<u64>(v);
}

std::vector<u32> series_divide(const std::vector<u32>& num, const std::vector<u32>& den,
                               const int degree) {
    // den[0] must be 1. Returns num/den modulo kMod up to x^degree.
    std::vector<u32> out(static_cast<std::size_t>(degree + 1), 0U);
    out[0] = num[0];
    for (int n = 1; n <= degree; ++n) {
        i128 acc = static_cast<i128>(num[static_cast<std::size_t>(n)]);
        for (int i = 1; i <= n; ++i) {
            acc -= static_cast<i128>(static_cast<u64>(den[static_cast<std::size_t>(i)]) *
                                     out[static_cast<std::size_t>(n - i)]);
        }
        out[static_cast<std::size_t>(n)] = static_cast<u32>(mod_norm(acc));
    }
    return out;
}

std::vector<u32> compute_F2_qseries(const int degree) {
    // F2(x) = 1/(1 - x^2/(1 - x^3/(1 - x^4/(...))))
    // Using the known identity (specializing the A047998 bivariate g.f. at y=x):
    //   F2(x) = (Sum_{n>=0} (-1)^n x^{n(n+2)} / (x;x)_n) / (Sum_{n>=0} (-1)^n x^{n(n+1)} / (x;x)_n),
    // where (x;x)_n = Prod_{k=1..n} (1 - x^k).
    std::vector<u32> den(static_cast<std::size_t>(degree + 1), 0U);
    std::vector<u32> num(static_cast<std::size_t>(degree + 1), 0U);

    std::vector<u32> invprod(static_cast<std::size_t>(degree + 1), 0U);
    invprod[0] = 1U;  // 1/(x;x)_0

    for (int n = 0;; ++n) {
        const int eD = n * (n + 1);
        const int eN = n * (n + 2);
        if (eD > degree && eN > degree) {
            break;
        }

        const bool neg = (n & 1) != 0;
        if (eD <= degree) {
            for (int i = 0; i + eD <= degree; ++i) {
                const u32 add = invprod[static_cast<std::size_t>(i)];
                u32& cell = den[static_cast<std::size_t>(i + eD)];
                if (!neg) {
                    const u64 s = static_cast<u64>(cell) + add;
                    cell = static_cast<u32>(s >= kMod ? s - kMod : s);
                } else {
                    cell = (cell >= add) ? static_cast<u32>(cell - add)
                                         : static_cast<u32>(cell + kMod - add);
                }
            }
        }
        if (eN <= degree) {
            for (int i = 0; i + eN <= degree; ++i) {
                const u32 add = invprod[static_cast<std::size_t>(i)];
                u32& cell = num[static_cast<std::size_t>(i + eN)];
                if (!neg) {
                    const u64 s = static_cast<u64>(cell) + add;
                    cell = static_cast<u32>(s >= kMod ? s - kMod : s);
                } else {
                    cell = (cell >= add) ? static_cast<u32>(cell - add)
                                         : static_cast<u32>(cell + kMod - add);
                }
            }
        }

        const int next = n + 1;
        if (next == 0) {
            continue;
        }
        for (int t = next; t <= degree; ++t) {
            const u64 s = static_cast<u64>(invprod[static_cast<std::size_t>(t)]) +
                          invprod[static_cast<std::size_t>(t - next)];
            invprod[static_cast<std::size_t>(t)] = static_cast<u32>(s % kMod);
        }
    }

    // den[0] == 1 so formal division is well-defined mod kMod.
    return series_divide(num, den, degree);
}

std::vector<u32> compute_F2_contfrac_small(const int degree) {
    // Direct series expansion of the continued fraction, for checkpointing only.
    const int k_max = degree + 5;
    std::vector<std::vector<u32>> F(static_cast<std::size_t>(k_max + 2));
    F[static_cast<std::size_t>(k_max + 1)] =
        std::vector<u32>(static_cast<std::size_t>(degree + 1), 0U);
    F[static_cast<std::size_t>(k_max + 1)][0] = 1U;

    for (int k = k_max; k >= 2; --k) {
        std::vector<u32> den(static_cast<std::size_t>(degree + 1), 0U);
        den[0] = 1U;
        const auto& next = F[static_cast<std::size_t>(k + 1)];
        for (int i = k; i <= degree; ++i) {
            const u32 v = next[static_cast<std::size_t>(i - k)];
            den[static_cast<std::size_t>(i)] = (v == 0U) ? 0U : static_cast<u32>(kMod - v);
        }

        // Invert den: inv[0]=1, inv[n] = -sum_{i=1..n} den[i]*inv[n-i]
        std::vector<u32> inv(static_cast<std::size_t>(degree + 1), 0U);
        inv[0] = 1U;
        for (int n = 1; n <= degree; ++n) {
            i128 acc = 0;
            for (int i = 1; i <= n; ++i) {
                acc -= static_cast<i128>(static_cast<u64>(den[static_cast<std::size_t>(i)]) *
                                         inv[static_cast<std::size_t>(n - i)]);
            }
            inv[static_cast<std::size_t>(n)] = static_cast<u32>(mod_norm(acc));
        }
        F[static_cast<std::size_t>(k)] = std::move(inv);
    }
    return F[2];
}

std::vector<u32> compute_t_from_F2(const std::vector<u32>& F2, const int degree) {
    // t(x) = x*(2*F2(x) - 1) / (1 - 2x(2*F2(x) - 1)) = T(x)/3
    std::vector<u32> A(static_cast<std::size_t>(degree + 1), 0U);     // A = 2F2 - 1
    std::vector<u32> den(static_cast<std::size_t>(degree + 1), 0U);   // 1 - 2xA
    std::vector<u32> num(static_cast<std::size_t>(degree + 1), 0U);   // xA

    for (int i = 0; i <= degree; ++i) {
        const u64 v = 2ULL * F2[static_cast<std::size_t>(i)];
        A[static_cast<std::size_t>(i)] = static_cast<u32>(v % kMod);
    }
    A[0] = (A[0] + static_cast<u32>(kMod - 1ULL)) % static_cast<u32>(kMod);

    den[0] = 1U;
    for (int i = 1; i <= degree; ++i) {
        const u64 v = 2ULL * A[static_cast<std::size_t>(i - 1)];
        den[static_cast<std::size_t>(i)] = static_cast<u32>(v == 0ULL ? 0ULL : (kMod - (v % kMod)) % kMod);
        num[static_cast<std::size_t>(i)] = A[static_cast<std::size_t>(i - 1)];
    }

    return series_divide(num, den, degree);
}

bool run_checkpoints() {
    constexpr int small_deg = 80;
    const auto f2_q = compute_F2_qseries(small_deg);
    const auto f2_cf = compute_F2_contfrac_small(small_deg);
    if (f2_q != f2_cf) {
        std::cerr << "Checkpoint failed: F2 q-series / contfrac mismatch\n";
        return false;
    }

    const auto t = compute_t_from_F2(f2_q, 20);
    const u64 T4 = (3ULL * t[4]) % kMod;
    const u64 T10 = (3ULL * t[10]) % kMod;
    if (T4 != 48ULL) {
        std::cerr << "Checkpoint failed: T(4)\n";
        return false;
    }
    if (T10 != 17'760ULL) {
        std::cerr << "Checkpoint failed: T(10)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr int n = 20'000;
    const auto f2 = compute_F2_qseries(n);
    const auto t = compute_t_from_F2(f2, n);
    const u64 answer = (3ULL * t[static_cast<std::size_t>(n)]) % kMod;
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}

