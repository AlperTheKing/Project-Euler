#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using i64 = long long;
using i128 = __int128_t;

static std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }
    bool neg = x < 0;
    if (neg) {
        x = -x;
    }
    std::string s;
    while (x > 0) {
        const int d = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    if (neg) {
        s.push_back('-');
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct MertensPrefix {
    int limit;
    std::vector<int> pref_mu;
    std::unordered_map<i64, i64> memo_mu;
    std::unordered_map<i64, i64> memo_odd_mu;

    explicit MertensPrefix(int lim) : limit(lim), pref_mu(lim + 1, 0) {
        std::vector<int> mu(lim + 1, 0);
        std::vector<int> primes;
        std::vector<bool> composite(lim + 1, false);

        mu[1] = 1;
        for (int i = 2; i <= lim; ++i) {
            if (!composite[i]) {
                primes.push_back(i);
                mu[i] = -1;
            }
            for (int p : primes) {
                const i64 v = 1LL * i * p;
                if (v > lim) {
                    break;
                }
                composite[static_cast<std::size_t>(v)] = true;
                if (i % p == 0) {
                    mu[static_cast<std::size_t>(v)] = 0;
                    break;
                }
                mu[static_cast<std::size_t>(v)] = -mu[i];
            }
        }

        for (int i = 1; i <= lim; ++i) {
            pref_mu[i] = pref_mu[i - 1] + mu[i];
        }
        memo_mu.reserve(1 << 16);
        memo_odd_mu.reserve(1 << 16);
    }

    i64 mu_prefix(i64 n) {
        if (n <= limit) {
            return pref_mu[static_cast<std::size_t>(n)];
        }
        const auto it = memo_mu.find(n);
        if (it != memo_mu.end()) {
            return it->second;
        }

        i64 ans = 1;
        i64 l = 2;
        while (l <= n) {
            const i64 q = n / l;
            const i64 r = n / q;
            ans -= (r - l + 1) * mu_prefix(q);
            l = r + 1;
        }

        memo_mu[n] = ans;
        return ans;
    }

    i64 odd_mu_prefix(i64 n) {
        if (n <= 0) {
            return 0;
        }
        const auto it = memo_odd_mu.find(n);
        if (it != memo_odd_mu.end()) {
            return it->second;
        }

        i64 ans = 0;
        i64 cur = n;
        while (cur > 0) {
            ans += mu_prefix(cur);
            cur >>= 1;
        }

        memo_odd_mu[n] = ans;
        return ans;
    }
};

static i128 sum_totients(i64 n, MertensPrefix& mp) {
    i128 s = 0;
    i64 l = 1;
    while (l <= n) {
        const i64 q = n / l;
        const i64 r = n / q;
        const i64 mu_seg = mp.mu_prefix(r) - mp.mu_prefix(l - 1);
        s += static_cast<i128>(mu_seg) * static_cast<i128>(q) * static_cast<i128>(q);
        l = r + 1;
    }
    return (s + 1) / 2;
}

static i64 f_odd_floor_sum(i64 m) {
    const i64 k = (m + 1) / 2;
    const i64 p = k / 2;
    if (k & 1LL) {
        return p * p;
    }
    return p * (p - 1);
}

static i128 solve_fast(i64 n) {
    const int sieve_limit = 2'000'000;
    MertensPrefix mp(sieve_limit);

    const i128 total = sum_totients(n, mp) - 1;

    i128 losing_half = 0;
    i64 l = 1;
    while (l <= n) {
        const i64 q = n / l;
        const i64 r = n / q;
        const i64 odd_mu_seg = mp.odd_mu_prefix(r) - mp.odd_mu_prefix(l - 1);
        losing_half += static_cast<i128>(odd_mu_seg) * static_cast<i128>(f_odd_floor_sum(q));
        l = r + 1;
    }

    const i128 losing = 2 * losing_half;
    return total - losing;
}

static bool has_losing_child(
    int a,
    int b,
    std::unordered_map<i64, bool>& memo,
    const std::function<bool(int, int)>& win
) {
    for (int u = 1; u < a; ++u) {
        const int m1 = b * u - 1;
        if (m1 % a == 0) {
            const int v = m1 / a;
            if (v > 0 && v < b) {
                if (!win(u, v)) {
                    return true;
                }
            }
        }

        const int m2 = b * u + 1;
        if (m2 % a == 0) {
            const int v = m2 / a;
            if (v > 0 && v < b) {
                if (!win(u, v)) {
                    return true;
                }
            }
        }
    }
    return false;
}

static i128 solve_bruteforce_game(int n) {
    std::unordered_map<i64, bool> memo;
    memo.reserve(1 << 16);

    std::function<bool(int, int)> win = [&](int a, int b) -> bool {
        const i64 key = (static_cast<i64>(a) << 32) | static_cast<unsigned int>(b);
        const auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        if (a == 1 || b == 1) {
            memo[key] = true;
            return true;
        }

        const bool w = has_losing_child(a, b, memo, win);
        memo[key] = w;
        return w;
    };

    i128 ans = 0;
    for (int a = 1; a < n; ++a) {
        for (int b = 1; a + b <= n; ++b) {
            if (std::gcd(a, b) != 1) {
                continue;
            }
            if (win(a, b)) {
                ++ans;
            }
        }
    }
    return ans;
}

int main() {
    assert(solve_fast(4) == 5);
    assert(solve_fast(100) == 2043);
    assert(solve_fast(40) == solve_bruteforce_game(40));

    const i128 ans = solve_fast(1'000'000'000LL);
    std::cout << to_string_i128(ans) << '\n';
    return 0;
}
