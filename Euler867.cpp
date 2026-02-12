#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;

static constexpr i64 MOD = 1'000'000'007LL;

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    while (e > 0) {
        if (e & 1LL) r = static_cast<i64>((static_cast<__int128>(r) * a) % MOD);
        a = static_cast<i64>((static_cast<__int128>(a) * a) % MOD);
        e >>= 1LL;
    }
    return r;
}

class Solver867 {
public:
    explicit Solver867(int n) : n_(n), max_len_(2 * n - 1), indep_(max_len_ + 1) {
        build_independent_masks();
    }

    i64 T(int n) {
        i64 ans = (2LL * R(n, n)) % MOD;
        if (n == 1) ans = (ans + MOD - 1) % MOD;
        return ans;
    }

private:
    int n_;
    int max_len_;
    std::vector<std::vector<int>> indep_;
    std::unordered_map<std::string, i64> count_cache_;
    std::unordered_map<int, i64> h_cache_;
    std::unordered_map<u64, i64> f_cache_;
    std::unordered_map<u64, i64> r_cache_;

    static u64 pair_key(int a, int b) {
        return (static_cast<u64>(static_cast<std::uint32_t>(a)) << 32) |
               static_cast<u64>(static_cast<std::uint32_t>(b));
    }

    void build_independent_masks() {
        for (int L = 0; L <= max_len_; ++L) {
            int lim = 1 << L;
            indep_[L].reserve(lim);
            for (int m = 0; m < lim; ++m) {
                if ((m & (m << 1)) == 0) indep_[L].push_back(m);
            }
        }
    }

    std::vector<i64> zeta_subset_sum(const std::vector<i64>& v, int L) const {
        std::vector<i64> out = v;
        int size = 1 << L;
        for (int bit = 0; bit < L; ++bit) {
            int step = 1 << bit;
            int block = step << 1;
            for (int start = 0; start < size; start += block) {
                int mid = start + step;
                int end = start + block;
                for (int m = mid; m < end; ++m) {
                    i64 x = out[m] + out[m - step];
                    out[m] = (x >= MOD ? x - MOD : x);
                }
            }
        }
        return out;
    }

    static std::string encode_lengths(const std::vector<int>& lens) {
        std::string key;
        key.reserve(lens.size() * 2 + 1);
        for (int x : lens) {
            key.push_back(static_cast<char>(x + 1));
            key.push_back('|');
        }
        return key;
    }

    i64 count_independent_sets(const std::vector<int>& row_lengths) {
        std::string key = encode_lengths(row_lengths);
        auto it = count_cache_.find(key);
        if (it != count_cache_.end()) return it->second;

        if (row_lengths.empty()) {
            count_cache_[key] = 1;
            return 1;
        }

        int L0 = row_lengths[0];
        std::vector<i64> dp(1 << L0, 0);
        for (int m : indep_[L0]) dp[m] = 1;

        for (std::size_t i = 1; i < row_lengths.size(); ++i) {
            int Lc = row_lengths[i - 1];
            int Ln = row_lengths[i];

            std::vector<i64> subs = zeta_subset_sum(dp, Lc);
            std::vector<i64> ndp(1 << Ln, 0);
            int full = (1 << Lc) - 1;

            if (Ln == Lc + 1) {
                for (int b : indep_[Ln]) {
                    int forb = (b | (b >> 1)) & full;
                    int allowed = full ^ forb;
                    ndp[b] = subs[allowed];
                }
            } else {
                for (int b : indep_[Ln]) {
                    int forb = (b | (b << 1)) & full;
                    int allowed = full ^ forb;
                    ndp[b] = subs[allowed];
                }
            }

            dp.swap(ndp);
        }

        int lastL = row_lengths.back();
        i64 total = 0;
        for (int m : indep_[lastL]) {
            total += dp[m];
            if (total >= MOD) total -= MOD;
        }

        count_cache_[key] = total;
        return total;
    }

    i64 H(int n) {
        auto it = h_cache_.find(n);
        if (it != h_cache_.end()) return it->second;

        std::vector<int> lens;
        lens.reserve(2 * n - 1);
        for (int x = n; x <= 2 * n - 1; ++x) lens.push_back(x);
        for (int x = 2 * n - 2; x >= n; --x) lens.push_back(x);

        i64 val = count_independent_sets(lens);
        h_cache_[n] = val;
        return val;
    }

    i64 F(int n, int h) {
        u64 key = pair_key(n, h);
        auto it = f_cache_.find(key);
        if (it != f_cache_.end()) return it->second;

        int rows = h - 1;
        std::vector<int> lens;
        lens.reserve(std::max(0, rows));
        for (int i = 0; i < rows; ++i) {
            lens.push_back(std::max(0, (n - 2) - i));
        }

        i64 val = count_independent_sets(lens);
        f_cache_[key] = val;
        return val;
    }

    i64 R(int u, int v) {
        u64 key = pair_key(u, v);
        auto it = r_cache_.find(key);
        if (it != r_cache_.end()) return it->second;

        i64 res;
        if (v == 0) {
            res = H(u);
        } else {
            res = (u == 1 && v == 1) ? 1 : 0;
            for (int w = 0; w < u; ++w) {
                i64 corner = F(u, u - w);
                i64 add = static_cast<i64>((static_cast<__int128>(R(v, w)) * mod_pow(corner, 6)) % MOD);
                res += add;
                if (res >= MOD) res -= MOD;
            }
        }

        r_cache_[key] = res;
        return res;
    }
};

int main() {
    Solver867 solver(10);

    assert(solver.T(1) == 5);
    assert(solver.T(2) == 48);

    std::cout << solver.T(10) << '\n';
    return 0;
}
