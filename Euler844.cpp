#include <algorithm>
#include <cassert>
#include <cstdint>
#include <deque>
#include <iostream>
#include <unordered_set>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kMod = 1'405'695'061ULL;

static inline u64 mod_add(u64 a, u64 b) {
    a += b;
    if (a >= kMod) a -= kMod;
    return a;
}

static inline u64 mod_sub(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

static inline u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

static u64 mod_pow(u64 a, u64 e) {
    u64 r = 1;
    while (e > 0) {
        if (e & 1ULL) r = mod_mul(r, a);
        a = mod_mul(a, a);
        e >>= 1ULL;
    }
    return r;
}

static u64 sum1_prefix(u64 n) {
    static const u64 inv2 = mod_pow(2, kMod - 2);
    const u64 a = n % kMod;
    const u64 b = (n + 1) % kMod;
    return mod_mul(mod_mul(a, b), inv2);
}

static u64 sum2_prefix(u64 n) {
    static const u64 inv6 = mod_pow(6, kMod - 2);
    const u64 a = n % kMod;
    const u64 b = (n + 1) % kMod;
    const u64 c = (2 * (n % kMod) + 1) % kMod;
    return mod_mul(mod_mul(mod_mul(a, b), c), inv6);
}

static u64 sum3_prefix(u64 n) {
    const u64 s1 = sum1_prefix(n);
    return mod_mul(s1, s1);
}

static u64 range_sum1(u64 l, u64 r) {
    if (l > r) return 0;
    return mod_sub(sum1_prefix(r), sum1_prefix(l - 1));
}

static u64 range_sum2(u64 l, u64 r) {
    if (l > r) return 0;
    return mod_sub(sum2_prefix(r), sum2_prefix(l - 1));
}

static u64 range_sum3(u64 l, u64 r) {
    if (l > r) return 0;
    return mod_sub(sum3_prefix(r), sum3_prefix(l - 1));
}

static u64 u2_value(u64 k) {
    return k * k - k - 1;
}

static u128 u3_value_u128(u64 k) {
    return static_cast<u128>(k) * k * k - static_cast<u128>(k) * k - 2 * static_cast<u128>(k) + 1;
}

static u128 p4_seed_u128(u64 k) {
    const u128 a = static_cast<u128>(k);
    return a * (k - 1) * (static_cast<u128>(k) * k - k - 1) - 1;
}

static u64 max_k_u2(u64 n) {
    u64 lo = 0, hi = 2'000'000'000ULL;
    while (lo < hi) {
        u64 mid = lo + (hi - lo + 1) / 2;
        if (u2_value(mid) <= n) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

static u64 max_k_u3(u64 n) {
    u64 lo = 0, hi = 2'000'000ULL;
    while (lo < hi) {
        u64 mid = lo + (hi - lo + 1) / 2;
        if (u3_value_u128(mid) <= static_cast<u128>(n)) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

static u64 max_k_seed4(u64 n) {
    u64 lo = 0, hi = 2'000'000ULL;
    while (lo < hi) {
        u64 mid = lo + (hi - lo + 1) / 2;
        if (p4_seed_u128(mid) <= static_cast<u128>(n)) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

struct VecHash {
    std::size_t operator()(const std::vector<u64>& v) const noexcept {
        std::size_t h = v.size();
        for (u64 x : v) {
            std::size_t y = std::hash<u64>{}(x);
            h ^= y + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

static bool capped_product_except(const std::vector<u64>& v, int skip, u64 cap, u64& out) {
    u128 prod = 1;
    for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        if (i == skip) continue;
        if (prod > static_cast<u128>(cap) / v[i]) return false;
        prod *= v[i];
    }
    out = static_cast<u64>(prod);
    return true;
}

static u64 exact_M_k(u64 k, u64 n) {
    std::deque<std::vector<u64>> q;
    std::unordered_set<std::vector<u64>, VecHash> seen;
    std::unordered_set<u64> numbers;

    q.push_back({});
    seen.insert({});
    numbers.insert(1);

    while (!q.empty()) {
        std::vector<u64> cur = q.front();
        q.pop_front();

        for (u64 x : cur) {
            if (x <= n) numbers.insert(x);
        }

        const u64 ones = k - static_cast<u64>(cur.size());

        if (ones > 0) {
            const u64 cap = (n + 1) / k;
            u64 prod = 1;
            if (capped_product_except(cur, -1, cap, prod)) {
                u64 y = k * prod - 1;
                if (y <= n) {
                    std::vector<u64> nxt = cur;
                    nxt.push_back(y);
                    auto it = std::upper_bound(nxt.begin(), nxt.end() - 1, y);
                    std::rotate(it, nxt.end() - 1, nxt.end());
                    if (seen.insert(nxt).second) q.push_back(std::move(nxt));
                }
            }
        }

        for (int i = 0; i < static_cast<int>(cur.size()); ++i) {
            if (i > 0 && cur[i] == cur[i - 1]) continue;
            u64 x = cur[i];
            const u64 cap = (n + x) / k;
            u64 prod_others = 1;
            if (!capped_product_except(cur, i, cap, prod_others)) continue;

            u64 y = k * prod_others - x;
            if (y == 0) continue;

            std::vector<u64> nxt;
            nxt.reserve(cur.size());
            for (int j = 0; j < static_cast<int>(cur.size()); ++j) {
                if (j != i) nxt.push_back(cur[j]);
            }
            if (y > 1) {
                auto it = std::upper_bound(nxt.begin(), nxt.end(), y);
                nxt.insert(it, y);
            }
            if (!nxt.empty() && nxt.back() > n) continue;
            if (seen.insert(nxt).second) q.push_back(std::move(nxt));
        }
    }

    u64 out = 0;
    for (u64 x : numbers) out = (out + (x % kMod)) % kMod;
    return out;
}

static u64 chain_sum_range(u64 l, u64 r, u64 n, u64 k2max, u64 k3max) {
    if (l > r) return 0;

    u64 ans = range_sum1(l, r);

    if (l <= k2max) {
        u64 rr = std::min(r, k2max);
        u64 cnt = (rr - l + 1) % kMod;
        u64 s2 = range_sum2(l, rr);
        u64 s1 = range_sum1(l, rr);
        u64 add = mod_sub(mod_sub(s2, s1), cnt);
        ans = mod_add(ans, add);
    }

    if (l <= k3max) {
        u64 rr = std::min(r, k3max);
        u64 cnt = (rr - l + 1) % kMod;
        u64 s3 = range_sum3(l, rr);
        u64 s2 = range_sum2(l, rr);
        u64 s1 = range_sum1(l, rr);
        u64 add = mod_sub(mod_sub(s3, s2), mod_mul(2, s1));
        add = mod_add(add, cnt);
        ans = mod_add(ans, add);
    }

    (void)n;
    return ans;
}

static u64 solve(u64 K, u64 N) {
    if (K < 3) return 0;

    const u64 k2max = max_k_u2(N);
    const u64 k3max = max_k_u3(N);
    const u64 k0 = max_k_seed4(N);

    const u64 exact_to = std::min(K, k0);
    u64 ans = 0;
    for (u64 k = 3; k <= exact_to; ++k) {
        ans = mod_add(ans, exact_M_k(k, N));
    }

    if (K > exact_to) {
        const u64 l = std::max<u64>(3, exact_to + 1);
        ans = mod_add(ans, chain_sum_range(l, K, N, k2max, k3max));
    }

    return ans;
}

int main() {
    assert(solve(4, 100) == 229U);
    assert(solve(10, 100'000'000ULL) == 2'383'369'980ULL % kMod);
    assert(exact_M_k(8, 100'000'000ULL) == 131'493'335ULL % kMod);

    std::cout << solve(1'000'000'000'000'000'000ULL, 1'000'000'000'000'000'000ULL) << '\n';
    return 0;
}
