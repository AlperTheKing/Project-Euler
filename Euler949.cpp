#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using i32 = int32_t;
using u64 = uint64_t;
using u128 = __uint128_t;

static constexpr u64 MOD = 1001001011ULL;

static inline u64 mod_add(u64 a, u64 b, u64 mod) {
    if (mod == 0) {
        return a + b;
    }
    a += b;
    if (a >= mod) a -= mod;
    return a;
}

static inline u64 mod_mul(u64 a, u64 b, u64 mod) {
    if (mod == 0) {
        return a * b;
    }
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

static inline i32 floor_div_pow2(i32 x, int s) {
    if (s == 0) return x;
    if (x >= 0) return static_cast<i32>(x >> s);
    const i32 ax = static_cast<i32>(-x);
    const i32 add = static_cast<i32>((1u << s) - 1u);
    return static_cast<i32>(-((ax + add) >> s));
}

static inline i32 ceil_div_pow2(i32 x, int s) {
    if (s == 0) return x;
    if (x >= 0) return static_cast<i32>((x + static_cast<i32>((1u << s) - 1u)) >> s);
    return static_cast<i32>(-((static_cast<i32>(-x)) >> s));
}

static i32 simplest_between(i32 u, i32 d, int e) {
    for (int m = 0; m <= e; ++m) {
        const int s = e - m;
        const i32 p_min = static_cast<i32>(floor_div_pow2(u, s) + 1);
        const i32 p_max = static_cast<i32>(ceil_div_pow2(d, s) - 1);
        if (p_min > p_max) continue;

        i32 p = 0;
        if (p_min > 0) {
            p = p_min;
        } else if (p_max < 0) {
            p = p_max;
        } else {
            p = 0;
        }

        if (m > 0 && p != 0 && (p & 1) == 0) {
            if (p + 1 <= p_max && ((p + 1) & 1) == 1) {
                ++p;
            } else if (p - 1 >= p_min && ((p - 1) & 1) == 1) {
                --p;
            }
        }

        return static_cast<i32>(static_cast<int64_t>(p) << s);
    }
    return 0;
}

struct UHotResult {
    std::vector<i32> u_full;
    std::vector<uint8_t> is_hot;
};

static UHotResult compute_u_hot(int n) {
    const int e = n;
    const i32 scale = static_cast<i32>(1u << e);
    const int total = static_cast<int>((1u << (n + 1)) - 1u);

    std::vector<i32> dp_u(total);
    std::vector<i32> dp_d(total);

    const int start1 = 1;
    dp_u[start1 + 0] = scale;
    dp_d[start1 + 0] = scale;
    dp_u[start1 + 1] = -scale;
    dp_d[start1 + 1] = -scale;

    std::vector<uint8_t> hot(static_cast<size_t>(1u << n), 0);

    for (int length = 2; length <= n; ++length) {
        const int size = 1u << length;
        const int start = (1u << length) - 1u;
        const bool final_layer = (length == n);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (size >= (1 << 15))
#endif
        for (int bits = 0; bits < size; ++bits) {
            i32 u_raw = std::numeric_limits<i32>::min();
            for (int s_len = 1; s_len < length; ++s_len) {
                const int suf = bits & ((1u << s_len) - 1u);
                const int idx = (1u << s_len) - 1u + suf;
                const i32 cand = dp_d[idx];
                if (cand > u_raw) u_raw = cand;
            }

            i32 d_raw = std::numeric_limits<i32>::max();
            for (int p_len = 1; p_len < length; ++p_len) {
                const int pre = bits >> (length - p_len);
                const int idx = (1u << p_len) - 1u + pre;
                const i32 cand = dp_u[idx];
                if (cand < d_raw) d_raw = cand;
            }

            const int idx = start + bits;
            if (u_raw < d_raw) {
                const i32 x = simplest_between(u_raw, d_raw, e);
                dp_u[idx] = x;
                dp_d[idx] = x;
                if (final_layer) hot[static_cast<size_t>(bits)] = 0;
            } else {
                dp_u[idx] = u_raw;
                dp_d[idx] = d_raw;
                if (final_layer) hot[static_cast<size_t>(bits)] = 1;
            }
        }
    }

    const int start_n = (1u << n) - 1u;
    const int full_size = 1u << n;

    UHotResult out;
    out.u_full.assign(dp_u.begin() + start_n, dp_u.begin() + start_n + full_size);
    out.is_hot = std::move(hot);
    return out;
}

using Dist = std::unordered_map<i32, u64>;

static Dist convolve(const Dist& a_in, const Dist& b_in, u64 mod) {
    if (a_in.empty() || b_in.empty()) return {};
    const Dist* a = &a_in;
    const Dist* b = &b_in;
    if (a->size() > b->size()) std::swap(a, b);

    Dist out;
    out.reserve(a->size() * std::min<size_t>(b->size(), 2048));

    for (const auto& [xa, ca] : *a) {
        for (const auto& [xb, cb] : *b) {
            const i32 key = static_cast<i32>(xa + xb);
            const u64 add = mod_mul(ca, cb, mod);
            auto it = out.find(key);
            if (it == out.end()) {
                out.emplace(key, add);
            } else {
                it->second = mod_add(it->second, add, mod);
            }
        }
    }
    return out;
}

static Dist pow_small(const Dist& hist, int t, u64 mod) {
    if (t == 0) return Dist{{0, 1}};
    Dist d = hist;
    for (int i = 1; i < t; ++i) {
        d = convolve(d, hist, mod);
    }
    return d;
}

static u64 count_sum_lt_zero(const Dist& a, const Dist& b, u64 mod) {
    std::vector<std::pair<i32, u64>> b_items;
    b_items.reserve(b.size());
    for (const auto& kv : b) b_items.push_back(kv);
    std::sort(b_items.begin(), b_items.end(),
              [](const auto& x, const auto& y) { return x.first < y.first; });

    std::vector<i32> b_sums;
    b_sums.reserve(b_items.size());
    for (const auto& [s, _] : b_items) b_sums.push_back(s);

    std::vector<u64> pref(b_items.size() + 1, 0);
    for (size_t i = 0; i < b_items.size(); ++i) {
        pref[i + 1] = mod_add(pref[i], b_items[i].second, mod);
    }

    u64 ans = 0;
    for (const auto& [sa, ca] : a) {
        const i32 target = static_cast<i32>(-sa);
        const auto it = std::lower_bound(b_sums.begin(), b_sums.end(), target);
        const size_t idx = static_cast<size_t>(it - b_sums.begin());
        ans = mod_add(ans, mod_mul(ca, pref[idx], mod), mod);
    }
    return ans;
}

static u64 count_sum_eq_zero(const Dist& a_in, const Dist& b_in, u64 mod) {
    const Dist* a = &a_in;
    const Dist* b = &b_in;
    if (a->size() > b->size()) std::swap(a, b);

    u64 ans = 0;
    for (const auto& [s, ca] : *a) {
        auto it = b->find(static_cast<i32>(-s));
        if (it == b->end()) continue;
        ans = mod_add(ans, mod_mul(ca, it->second, mod), mod);
    }
    return ans;
}

static u64 G(int n, int k, u64 mod) {
    if (k <= 0 || (k % 2) == 0) throw std::invalid_argument("k must be positive and odd");
    if (n <= 0) throw std::invalid_argument("n must be positive");

    const auto [u_full, is_hot] = compute_u_hot(n);

    Dist hist_all;
    Dist hist_cold;
    hist_all.reserve(1u << 12);
    hist_cold.reserve(1u << 12);

    for (size_t i = 0; i < u_full.size(); ++i) {
        const i32 v = u_full[i];
        {
            auto it = hist_all.find(v);
            if (it == hist_all.end()) {
                hist_all.emplace(v, (mod == 0) ? 1ULL : 1ULL % mod);
            } else {
                it->second = mod_add(it->second, 1, mod);
            }
        }
        if (is_hot[i] == 0) {
            auto it = hist_cold.find(v);
            if (it == hist_cold.end()) {
                hist_cold.emplace(v, (mod == 0) ? 1ULL : 1ULL % mod);
            } else {
                it->second = mod_add(it->second, 1, mod);
            }
        }
    }

    const int a = k / 2;
    const int b = k - a;

    const Dist dist_a = pow_small(hist_all, a, mod);
    const Dist dist_b = pow_small(hist_all, b, mod);
    const u64 neg = count_sum_lt_zero(dist_a, dist_b, mod);

    const Dist cold_a = pow_small(hist_cold, a, mod);
    const Dist cold_b = pow_small(hist_cold, b, mod);
    const u64 zero_cold = count_sum_eq_zero(cold_a, cold_b, mod);

    return mod_add(neg, zero_cold, mod);
}

int main() {
#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << "\n";
#endif

    const u64 g23 = G(2, 3, 0);
    std::cout << "G(2,3) = " << g23 << " (expected 14)\n";
    if (g23 != 14) return 1;

    const u64 g43 = G(4, 3, 0);
    std::cout << "G(4,3) = " << g43 << " (expected 496)\n";
    if (g43 != 496) return 1;

    const u64 g85 = G(8, 5, 0);
    std::cout << "G(8,5) = " << g85 << " (expected 26359197010)\n";
    if (g85 != 26359197010ULL) return 1;

    const u64 ans = G(20, 7, MOD);
    std::cout << "G(20,7) mod " << MOD << " = " << ans << "\n";
    return 0;
}
