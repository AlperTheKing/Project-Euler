#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using i64 = long long;
using u64 = unsigned long long;

constexpr i64 kMod = 1000000007LL;
constexpr i64 kDefaultN = 100000000000000LL;

i64 mod_pow(i64 base, i64 exp) {
    i64 res = 1 % kMod;
    i64 cur = base % kMod;
    while (exp > 0) {
        if (exp & 1LL) res = (res * cur) % kMod;
        cur = (cur * cur) % kMod;
        exp >>= 1LL;
    }
    return res;
}

struct Precomputed {
    int max_deg = 0;
    std::vector<std::vector<i64>> binom;
    std::vector<std::vector<i64>> sum_pow;
};

std::vector<int> to_digits(i64 n, int base) {
    if (n == 0) return {0};
    std::vector<int> digits;
    while (n > 0) {
        digits.push_back(static_cast<int>(n % base));
        n /= base;
    }
    return digits;
}

Precomputed build_precomputed(int max_deg) {
    Precomputed pre;
    pre.max_deg = max_deg;
    const int size = max_deg + 2;

    std::vector<std::vector<i64>> binom(size, std::vector<i64>(size, 0));
    for (int n = 0; n < size; ++n) {
        binom[n][0] = 1;
        binom[n][n] = 1;
        for (int k = 1; k < n; ++k) {
            binom[n][k] = binom[n - 1][k - 1] + binom[n - 1][k];
            if (binom[n][k] >= kMod) binom[n][k] -= kMod;
        }
    }

    std::vector<std::vector<i64>> stirling(size, std::vector<i64>(size, 0));
    stirling[0][0] = 1;
    for (int n = 1; n < size; ++n) {
        for (int k = 1; k <= n; ++k) {
            stirling[n][k] = (stirling[n - 1][k - 1] + k * stirling[n - 1][k]) % kMod;
        }
    }

    std::vector<i64> fact(size + 1, 1);
    for (int i = 1; i < static_cast<int>(fact.size()); ++i) {
        fact[i] = (fact[i - 1] * i) % kMod;
    }

    std::vector<std::vector<i64>> poly_binom(size);
    poly_binom[0] = {1};
    for (int m = 1; m < size; ++m) {
        const auto& prev = poly_binom[m - 1];
        std::vector<i64> cur(prev.size() + 1, 0);
        for (int p = 0; p < static_cast<int>(prev.size()); ++p) {
            const i64 a = prev[p];
            const i64 dec = (a * (m - 1)) % kMod;
            cur[p] -= dec;
            if (cur[p] < 0) cur[p] += kMod;
            cur[p + 1] += a;
            if (cur[p + 1] >= kMod) cur[p + 1] -= kMod;
        }
        const i64 inv_m = mod_pow(m, kMod - 2);
        for (auto& v : cur) {
            v = (v * inv_m) % kMod;
        }
        poly_binom[m] = std::move(cur);
    }

    std::vector<std::vector<i64>> poly_binom_shift(size);
    for (int m = 0; m < size; ++m) {
        const auto& poly = poly_binom[m];
        std::vector<i64> shift(poly.size(), 0);
        for (int p = 0; p < static_cast<int>(poly.size()); ++p) {
            const i64 a = poly[p];
            if (a == 0) continue;
            for (int q = 0; q <= p; ++q) {
                i64 add = (a * binom[p][q]) % kMod;
                shift[q] += add;
                if (shift[q] >= kMod) shift[q] -= kMod;
            }
        }
        poly_binom_shift[m] = std::move(shift);
    }

    std::vector<std::vector<i64>> sum_pow(max_deg + 1);
    for (int p = 0; p <= max_deg; ++p) {
        std::vector<i64> poly(p + 2, 0);
        for (int j = 0; j <= p; ++j) {
            if (stirling[p][j] == 0) continue;
            const i64 coeff = (stirling[p][j] * fact[j]) % kMod;
            const auto& base = poly_binom_shift[j + 1];  // C(x+1, j+1)
            for (int idx = 0; idx < static_cast<int>(base.size()); ++idx) {
                i64 add = (coeff * base[idx]) % kMod;
                poly[idx] += add;
                if (poly[idx] >= kMod) poly[idx] -= kMod;
            }
        }
        sum_pow[p] = std::move(poly);
    }

    pre.binom = std::move(binom);
    pre.sum_pow = std::move(sum_pow);
    return pre;
}

std::vector<i64> prefix_sum_poly(const std::vector<i64>& poly, const Precomputed& pre) {
    std::vector<i64> res(poly.size() + 1, 0);
    for (int p = 0; p < static_cast<int>(poly.size()); ++p) {
        const i64 a = poly[p];
        if (a == 0) continue;
        const auto& base = pre.sum_pow[p];
        for (int i = 0; i < static_cast<int>(base.size()); ++i) {
            i64 add = (a * base[i]) % kMod;
            res[i] += add;
            if (res[i] >= kMod) res[i] -= kMod;
        }
    }
    return res;
}

std::vector<i64> substitute_poly(const std::vector<i64>& poly, int k, int s,
                                 const Precomputed& pre, const std::vector<i64>& pow_k) {
    const int d = static_cast<int>(poly.size()) - 1;
    std::vector<i64> pow_s(d + 1, 1);
    for (int i = 1; i <= d; ++i) {
        pow_s[i] = (pow_s[i - 1] * s) % kMod;
    }

    std::vector<i64> res(d + 1, 0);
    for (int p = 0; p <= d; ++p) {
        const i64 a = poly[p];
        if (a == 0) continue;
        for (int q = 0; q <= p; ++q) {
            i64 term = (a * pre.binom[p][q]) % kMod;
            term = (term * pow_k[q]) % kMod;
            term = (term * pow_s[p - q]) % kMod;
            res[q] += term;
            if (res[q] >= kMod) res[q] -= kMod;
        }
    }
    return res;
}

i64 compute_f(i64 n, int k, const Precomputed& pre) {
    if (n == 0) return 1;
    // Carry-based digit DP: polynomials encode counts as a function of the next carry.
    std::vector<int> digits = to_digits(n, k);
    const int L = static_cast<int>(digits.size()) - 1;

    std::vector<i64> pow_k(pre.max_deg + 2, 1);
    for (int i = 1; i < static_cast<int>(pow_k.size()); ++i) {
        pow_k[i] = (pow_k[i - 1] * k) % kMod;
    }

    std::vector<i64> free_poly(1, k % kMod);
    std::vector<i64> tight_poly(1, (digits[0] + 1) % kMod);
    if (L == 0) return tight_poly[0];

    for (int j = 1; j <= L; ++j) {
        std::vector<i64> sum_free = prefix_sum_poly(free_poly, pre);
        std::vector<i64> sum_tight = prefix_sum_poly(tight_poly, pre);

        std::vector<i64> new_free(sum_free.size(), 0);
        for (int s = 0; s < k; ++s) {
            std::vector<i64> sub = substitute_poly(sum_free, k, s, pre, pow_k);
            for (int i = 0; i < static_cast<int>(sub.size()); ++i) {
                new_free[i] += sub[i];
                if (new_free[i] >= kMod) new_free[i] -= kMod;
            }
        }

        std::vector<i64> new_tight(sum_free.size(), 0);
        const int dj = digits[j];
        for (int s = 0; s < dj; ++s) {
            std::vector<i64> sub = substitute_poly(sum_free, k, s, pre, pow_k);
            for (int i = 0; i < static_cast<int>(sub.size()); ++i) {
                new_tight[i] += sub[i];
                if (new_tight[i] >= kMod) new_tight[i] -= kMod;
            }
        }
        std::vector<i64> sub_tight = substitute_poly(sum_tight, k, dj, pre, pow_k);
        for (int i = 0; i < static_cast<int>(sub_tight.size()); ++i) {
            new_tight[i] += sub_tight[i];
            if (new_tight[i] >= kMod) new_tight[i] -= kMod;
        }

        free_poly = std::move(new_free);
        tight_poly = std::move(new_tight);
    }

    return tight_poly[0] % kMod;
}

u64 naive_f(int k, int n) {
    std::vector<u64> g(n + 1, 0);
    g[0] = 1;
    for (int i = 1; i <= n; ++i) {
        g[i] = g[i - 1] + g[i / k];
    }
    return g[n];
}

void run_validation(const Precomputed& pre) {
    struct Example { int k; int n; u64 expected; };
    const Example examples[] = {
        {5, 10, 18ULL},
        {7, 100, 1003ULL},
        {2, 1000, 264830889564ULL},
    };
    for (const auto& ex : examples) {
        const i64 got = compute_f(ex.n, ex.k, pre);
        const i64 want = static_cast<i64>(ex.expected % kMod);
        if (got != want) {
            std::cerr << "Validation failed for k=" << ex.k << " n=" << ex.n
                      << ": got " << got << " expected " << want << '\n';
            std::exit(1);
        }
    }

    for (int k = 2; k <= 10; ++k) {
        const int n = 200;
        const u64 exact = naive_f(k, n);
        const i64 got = compute_f(n, k, pre);
        if (got != static_cast<i64>(exact % kMod)) {
            std::cerr << "Validation failed for k=" << k << " n=" << n
                      << ": got " << got << " expected " << (exact % kMod) << '\n';
            std::exit(1);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    i64 n = kDefaultN;
    if (argc > 1) {
        n = std::strtoll(argv[1], nullptr, 10);
    }

    const i64 max_n = std::max<i64>(n, 1000);
    const int max_deg = static_cast<int>(to_digits(max_n, 2).size()) + 2;
    const Precomputed pre = build_precomputed(max_deg);

    run_validation(pre);

    std::vector<std::future<i64>> futures;
    futures.reserve(9);
    for (int k = 2; k <= 10; ++k) {
        futures.emplace_back(std::async(std::launch::async, [&, k]() {
            return compute_f(n, k, pre);
        }));
    }

    i64 total = 0;
    for (auto& fut : futures) {
        i64 val = fut.get();
        total += val;
        if (total >= kMod) total -= kMod;
    }

    std::cout << total % kMod << '\n';
    return 0;
}
