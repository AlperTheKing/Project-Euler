#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;

static constexpr u64 MOD = 987654321ULL;

static u64 mod_pow(u64 a, u64 e) {
    a %= MOD;
    u64 r = 1 % MOD;
    while (e) {
        if (e & 1ULL) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1ULL;
    }
    return r;
}

static i64 ext_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0, y1 = 0;
    i64 g = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static u64 mod_inv(u64 a) {
    i64 x = 0, y = 0;
    i64 g = ext_gcd(static_cast<i64>(a % MOD), static_cast<i64>(MOD), x, y);
    assert(g == 1);
    i64 res = x % static_cast<i64>(MOD);
    if (res < 0) res += static_cast<i64>(MOD);
    return static_cast<u64>(res);
}

static void fwht_xor(std::vector<u64>& a, bool inverse) {
    const std::size_t n = a.size();
    for (std::size_t len = 1; len < n; len <<= 1U) {
        for (std::size_t i = 0; i < n; i += (len << 1U)) {
            for (std::size_t j = 0; j < len; ++j) {
                u64 u = a[i + j];
                u64 v = a[i + j + len];
                u64 x = u + v;
                if (x >= MOD) x -= MOD;
                u64 y = (u >= v) ? (u - v) : (u + MOD - v);
                a[i + j] = x;
                a[i + j + len] = y;
            }
        }
    }
    if (inverse) {
        const u64 inv_n = mod_inv(static_cast<u64>(n));
        for (u64& x : a) x = (x * inv_n) % MOD;
    }
}

static u64 encode_pattern(const int* exps, int len) {
    // Encode non-increasing exponent multiset into 64-bit:
    // top 4 bits: length (<=15), then 6 bits per exponent (<=63).
    u64 key = (static_cast<u64>(len) << 60U);
    for (int i = 0; i < len; ++i) key |= (static_cast<u64>(exps[i]) << (6U * i));
    return key;
}

static void decode_pattern(u64 key, int* exps, int& len) {
    len = static_cast<int>(key >> 60U);
    for (int i = 0; i < len; ++i) exps[i] = static_cast<int>((key >> (6U * i)) & 63ULL);
}

struct GrundyComputer {
    std::unordered_map<u64, u32> memo;

    GrundyComputer() { memo.reserve(1024); }

    u32 grundy(u64 key) {
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;

        int a[12];
        int len = 0;
        decode_pattern(key, a, len);

        // Enumerate all proper divisors d of n (excluding 1 and n itself) in exponent space.
        int cur[12];
        std::vector<u32> div_nims;
        div_nims.reserve(256);

        std::function<void(int, bool, bool)> dfs = [&](int idx, bool any_pos, bool all_eq) {
            if (idx == len) {
                if (!any_pos) return;        // d=1
                if (all_eq) return;          // d=n
                int b[12];
                int blen = 0;
                for (int i = 0; i < len; ++i) {
                    if (cur[i] > 0) b[blen++] = cur[i];
                }
                // insertion sort descending, blen is small
                for (int i = 1; i < blen; ++i) {
                    int v = b[i];
                    int j = i;
                    while (j > 0 && b[j - 1] < v) {
                        b[j] = b[j - 1];
                        --j;
                    }
                    b[j] = v;
                }
                const u64 dkey = encode_pattern(b, blen);
                div_nims.push_back(grundy(dkey));
                return;
            }
            for (int e = 0; e <= a[idx]; ++e) {
                cur[idx] = e;
                dfs(idx + 1, any_pos || (e > 0), all_eq && (e == a[idx]));
            }
        };
        dfs(0, false, true);

        if (div_nims.empty()) {
            memo.emplace(key, 0U);
            return 0U;
        }

        std::sort(div_nims.begin(), div_nims.end());
        div_nims.erase(std::unique(div_nims.begin(), div_nims.end()), div_nims.end());

        const u32 mx = div_nims.back();
        std::vector<std::uint8_t> present(static_cast<std::size_t>(mx + 1), 0);
        for (u32 v : div_nims) present[static_cast<std::size_t>(v)] = 1;

        u32 mex = 0;
        while (true) {
            bool ok = false;
            for (u32 v : div_nims) {
                const u32 t = v ^ mex;
                if (t <= mx && present[static_cast<std::size_t>(t)]) {
                    ok = true;
                    break;
                }
            }
            if (!ok) break;
            ++mex;
        }

        memo.emplace(key, mex);
        return mex;
    }
};

static std::vector<u32> build_spf(u32 n) {
    std::vector<u32> spf(static_cast<std::size_t>(n + 1), 0U);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    for (u32 i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (u32 p : primes) {
            const u64 v = static_cast<u64>(p) * static_cast<u64>(i);
            if (v > n) break;
            spf[static_cast<std::size_t>(v)] = p;
            if (p == spf[i]) break;
        }
    }
    return spf;
}

static u64 winning_positions(u32 n, u64 k, GrundyComputer& gc) {
    const std::vector<u32> spf = build_spf(n);

    std::vector<u64> counts(1, 0ULL);
    u32 max_g = 0;

    for (u32 x = 2; x <= n; ++x) {
        u32 t = x;
        int exps[12];
        int len = 0;
        while (t > 1) {
            const u32 p = spf[t];
            int c = 0;
            while (t > 1 && spf[t] == p) {
                t /= p;
                ++c;
            }
            exps[len++] = c;
        }
        // insertion sort descending
        for (int i = 1; i < len; ++i) {
            int v = exps[i];
            int j = i;
            while (j > 0 && exps[j - 1] < v) {
                exps[j] = exps[j - 1];
                --j;
            }
            exps[j] = v;
        }
        const u64 key = encode_pattern(exps, len);
        const u32 g = gc.grundy(key);
        if (g >= counts.size()) counts.resize(static_cast<std::size_t>(g + 1), 0ULL);
        ++counts[static_cast<std::size_t>(g)];
        if (g > max_g) max_g = g;
    }

    std::size_t sz = 1;
    while (sz <= static_cast<std::size_t>(max_g)) sz <<= 1U;
    std::vector<u64> a(sz, 0ULL);
    for (std::size_t i = 0; i < counts.size(); ++i) a[i] = counts[i] % MOD;

    fwht_xor(a, false);
    for (u64& x : a) x = mod_pow(x, k);
    fwht_xor(a, true);

    const u64 losing = a[0] % MOD;
    const u64 total = mod_pow(static_cast<u64>(n - 1U), k);
    return (total + MOD - losing) % MOD;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    GrundyComputer gc;

    // Validation point from the problem statement.
    assert(winning_positions(10, 5, gc) == 40085ULL);

    std::cout << winning_positions(10'000'000U, 1'000'000'000'000ULL, gc) % MOD << '\n';
    return 0;
}
