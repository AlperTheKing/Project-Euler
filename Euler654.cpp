#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_add(u64 a, u64 b) {
    a += b;
    if (a >= kMod) a -= kMod;
    return a;
}

u64 mod_sub(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

u64 mod_mul(u64 a, u64 b) {
    return (a * b) % kMod;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) result = mod_mul(result, base);
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

std::vector<u64> build_sequence(int n, int terms) {
    const int states = n - 1;
    std::vector<u64> vec(static_cast<std::size_t>(states), 1);
    std::vector<u64> prefix(static_cast<std::size_t>(states), 0);
    std::vector<u64> next_vec(static_cast<std::size_t>(states), 0);

    std::vector<u64> seq;
    seq.reserve(static_cast<std::size_t>(terms));

    for (int step = 1; step <= terms; ++step) {
        u64 total = 0;
        for (u64 x : vec) {
            total += x;
            if (total >= kMod) total -= kMod;
        }
        seq.push_back(total);

        u64 run = 0;
        for (int i = 0; i < states; ++i) {
            run += vec[static_cast<std::size_t>(i)];
            if (run >= kMod) run -= kMod;
            prefix[static_cast<std::size_t>(i)] = run;
        }

        for (int j = 0; j < states; ++j) {
            next_vec[static_cast<std::size_t>(j)] =
                prefix[static_cast<std::size_t>(states - 1 - j)];
        }
        vec.swap(next_vec);
    }

    return seq;
}

std::vector<u64> berlekamp_massey(const std::vector<u64>& s) {
    std::vector<u64> C{1}, B{1};
    int L = 0;
    int m = 1;
    u64 b = 1;

    for (int n = 0; n < static_cast<int>(s.size()); ++n) {
        u64 d = s[static_cast<std::size_t>(n)];
        for (int i = 1; i <= L; ++i) {
            d = (d + mod_mul(C[static_cast<std::size_t>(i)], s[static_cast<std::size_t>(n - i)])) % kMod;
        }

        if (d == 0) {
            ++m;
            continue;
        }

        const std::vector<u64> T = C;
        const u64 coef = mod_mul(d, mod_pow(b, kMod - 2));

        if (C.size() < B.size() + static_cast<std::size_t>(m)) {
            C.resize(B.size() + static_cast<std::size_t>(m), 0);
        }

        for (std::size_t i = 0; i < B.size(); ++i) {
            const std::size_t idx = i + static_cast<std::size_t>(m);
            C[idx] = mod_sub(C[idx], mod_mul(coef, B[i]));
        }

        if (2 * L <= n) {
            L = n + 1 - L;
            B = T;
            b = d;
            m = 1;
        } else {
            ++m;
        }
    }

    std::vector<u64> rec(static_cast<std::size_t>(L), 0);
    for (int i = 1; i <= L; ++i) {
        rec[static_cast<std::size_t>(i - 1)] = mod_sub(0, C[static_cast<std::size_t>(i)]);
    }
    return rec;
}

std::vector<u64> combine_poly(const std::vector<u64>& a,
                              const std::vector<u64>& b,
                              const std::vector<u64>& rec) {
    const int k = static_cast<int>(rec.size());
    std::vector<u64> tmp(static_cast<std::size_t>(2 * k - 1), 0);

    for (int i = 0; i < k; ++i) {
        if (a[static_cast<std::size_t>(i)] == 0) continue;
        for (int j = 0; j < k; ++j) {
            if (b[static_cast<std::size_t>(j)] == 0) continue;
            const u64 add =
                mod_mul(a[static_cast<std::size_t>(i)], b[static_cast<std::size_t>(j)]);
            u64& cell = tmp[static_cast<std::size_t>(i + j)];
            cell += add;
            if (cell >= kMod) cell -= kMod;
        }
    }

    for (int i = 2 * k - 2; i >= k; --i) {
        const u64 x = tmp[static_cast<std::size_t>(i)];
        if (x == 0) continue;
        for (int t = 1; t <= k; ++t) {
            const int idx = i - t;
            const u64 add = mod_mul(x, rec[static_cast<std::size_t>(t - 1)]);
            u64& cell = tmp[static_cast<std::size_t>(idx)];
            cell += add;
            if (cell >= kMod) cell -= kMod;
        }
    }

    tmp.resize(static_cast<std::size_t>(k));
    return tmp;
}

u64 linear_recurrence_nth(const std::vector<u64>& init,
                          const std::vector<u64>& rec,
                          u64 n_one_based) {
    const int k = static_cast<int>(rec.size());
    if (n_one_based <= static_cast<u64>(k)) return init[static_cast<std::size_t>(n_one_based - 1)];

    u64 exp = n_one_based - 1;
    std::vector<u64> result(static_cast<std::size_t>(k), 0);
    result[0] = 1;

    std::vector<u64> base(static_cast<std::size_t>(k), 0);
    if (k == 1) {
        base[0] = rec[0];
    } else {
        base[1] = 1;
    }

    while (exp > 0) {
        if (exp & 1ULL) result = combine_poly(result, base, rec);
        exp >>= 1ULL;
        if (exp > 0) base = combine_poly(base, base, rec);
    }

    u64 ans = 0;
    for (int i = 0; i < k; ++i) {
        ans = (ans + mod_mul(result[static_cast<std::size_t>(i)], init[static_cast<std::size_t>(i)])) % kMod;
    }
    return ans;
}

u64 solve_case(int n, u64 m) {
    const int need_terms = 2 * (n - 1) + 8;
    const std::vector<u64> seq = build_sequence(n, need_terms);
    if (m <= seq.size()) return seq[static_cast<std::size_t>(m - 1)];

    const std::vector<u64> rec = berlekamp_massey(seq);
    const int k = static_cast<int>(rec.size());
    std::vector<u64> init(seq.begin(), seq.begin() + k);
    return linear_recurrence_nth(init, rec, m);
}

}  // namespace

int main() {
    assert(solve_case(3, 4) == 8);
    assert(solve_case(5, 5) == 246);
    assert(solve_case(10, 100) == 862820094);
    assert(solve_case(100, 10) == 782136797);

    std::cout << solve_case(5000, 1'000'000'000'000ULL) << "\n";
    return 0;
}
