#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;

static u64 pow_limited(u64 base, int exp, u64 limit) {
    u64 res = 1;
    for (int i = 0; i < exp; ++i) {
        if (res > limit / base) return limit + 1;
        res *= base;
    }
    return res;
}

static u64 iroot(u64 n, int k) {
    if (k == 1 || n <= 1) return n;
    long double x = std::powl(static_cast<long double>(n), 1.0L / static_cast<long double>(k));
    u64 r = static_cast<u64>(x);
    while (pow_limited(r + 1, k, n) <= n) ++r;
    while (pow_limited(r, k, n) > n) --r;
    return r;
}

class PrimePi {
public:
    PrimePi(int limit = 1'000'000) : limit_(limit), pi_small_(limit + 1, 0) {
        std::vector<bool> is_prime(limit_ + 1, true);
        is_prime[0] = is_prime[1] = false;
        for (int p = 2; p * p <= limit_; ++p) {
            if (!is_prime[p]) continue;
            for (int q = p * p; q <= limit_; q += p) is_prime[q] = false;
        }
        for (int i = 2; i <= limit_; ++i) {
            if (is_prime[i]) primes_.push_back(i);
            pi_small_[i] = pi_small_[i - 1] + (is_prime[i] ? 1 : 0);
        }
    }

    const std::vector<int>& primes() const { return primes_; }

    u64 pi(u64 x) {
        if (x <= static_cast<u64>(limit_)) return static_cast<u64>(pi_small_[static_cast<std::size_t>(x)]);
        auto it = pi_cache_.find(x);
        if (it != pi_cache_.end()) return it->second;

        u64 a = pi(iroot(x, 4));
        u64 b = pi(iroot(x, 2));
        u64 c = pi(iroot(x, 3));

        u64 sum = phi(x, static_cast<int>(a)) + (b + a - 2) * (b - a + 1) / 2;

        for (u64 i = a + 1; i <= b; ++i) {
            u64 w = x / static_cast<u64>(primes_[static_cast<std::size_t>(i - 1)]);
            sum -= pi(w);
            if (i <= c) {
                u64 lim = pi(iroot(w, 2));
                for (u64 j = i; j <= lim; ++j) {
                    u64 pw = static_cast<u64>(primes_[static_cast<std::size_t>(j - 1)]);
                    sum -= pi(w / pw) - (j - 1);
                }
            }
        }

        pi_cache_[x] = sum;
        return sum;
    }

private:
    int limit_;
    std::vector<int> primes_;
    std::vector<int> pi_small_;
    std::unordered_map<u64, u64> phi_cache_;
    std::unordered_map<u64, u64> pi_cache_;

    static u64 phi_key(u64 x, int s) {
        return (x << 9) ^ static_cast<u64>(s);
    }

    u64 phi(u64 x, int s) {
        if (s == 0) return x;
        if (s == 1) return x - x / 2;

        u64 key = phi_key(x, s);
        auto it = phi_cache_.find(key);
        if (it != phi_cache_.end()) return it->second;

        u64 res = phi(x, s - 1) - phi(x / static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]), s - 1);
        phi_cache_[key] = res;
        return res;
    }
};

class Solver {
public:
    explicit Solver(u64 N) : N_(N), prime_pi_(1'000'000), primes_(prime_pi_.primes()) {}

    u64 Q(int k) {
        int D = 2 * k;
        auto seqs = exponent_sequences(D);
        u64 total = 0;
        for (const auto& seq : seqs) total += count_sequence(seq);
        return total;
    }

private:
    u64 N_;
    PrimePi prime_pi_;
    const std::vector<int>& primes_;

    static int c_value(int a) {
        return (a & 1) ? (a + 1) : a;
    }

    std::vector<std::vector<int>> exponent_sequences(int D) {
        static const int cvals[10] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
        std::vector<std::vector<int>> out;
        std::vector<int> cur;

        std::function<void(int)> dfs = [&](int rem) {
            if (rem == 1) {
                if (!cur.empty()) out.push_back(cur);
                return;
            }
            for (int c : cvals) {
                if (rem % c != 0) continue;
                cur.push_back(c - 1);
                dfs(rem / c);
                cur.pop_back();
                cur.push_back(c);
                dfs(rem / c);
                cur.pop_back();
            }
        };

        dfs(D);
        return out;
    }

    u64 count_sequence(const std::vector<int>& exps) {
        const int r = static_cast<int>(exps.size());
        std::vector<int> suffix(r + 1, 0);
        for (int i = r - 1; i >= 0; --i) suffix[i] = suffix[i + 1] + exps[i];

        std::unordered_map<u64, u64> memo;
        memo.reserve(1 << 14);

        std::function<u64(int, int, u64)> dfs = [&](int pos, int start_idx, u64 limit) -> u64 {
            if (pos == r) return 1;
            u64 key = (static_cast<u64>(pos) << 61) ^ (static_cast<u64>(start_idx) << 44) ^ limit;
            auto it = memo.find(key);
            if (it != memo.end()) return it->second;

            u64 ans = 0;
            if (pos == r - 1) {
                u64 mx = iroot(limit, exps[pos]);
                u64 cnt = prime_pi_.pi(mx);
                ans = (cnt > static_cast<u64>(start_idx)) ? (cnt - static_cast<u64>(start_idx)) : 0;
            } else {
                u64 mx = iroot(limit, suffix[pos]);
                auto it_hi = std::upper_bound(primes_.begin() + start_idx, primes_.end(), static_cast<int>(mx));
                int hi = static_cast<int>(it_hi - primes_.begin());

                for (int i = start_idx; i < hi; ++i) {
                    u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(i)]);
                    u64 pw = pow_limited(p, exps[pos], limit);
                    if (pw > limit) break;
                    ans += dfs(pos + 1, i + 1, limit / pw);
                }
            }

            memo[key] = ans;
            return ans;
        };

        return dfs(0, 0, N_);
    }
};

int main() {
    {
        Solver s(100);
        assert(s.Q(2) == 51);
    }
    {
        Solver s(1'000'000);
        assert(s.Q(6) == 6'189);
    }

    Solver s(1'000'000'000'000ULL);
    u64 ans = 0;
    for (int k = 2; k <= 10; ++k) ans += s.Q(k);
    std::cout << ans << '\n';
    return 0;
}
