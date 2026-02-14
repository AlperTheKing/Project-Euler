#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) ++r;
    while (r > x / r) --r;
    return r;
}

std::vector<int> sieve_primes(int limit) {
    if (limit < 2) return {};
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = is_prime[1] = 0U;

    for (int i = 2; 1LL * i * i <= limit; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) continue;
        for (int j = i * i; j <= limit; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0U;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::max(1.0, std::log(static_cast<double>(limit)))));
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string out;
    while (value > 0) {
        out.push_back(static_cast<char>('0' + static_cast<int>(value % 10)));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

class Solver708 {
  public:
    explicit Solver708(u64 n) : n_(n), primes_(sieve_primes(static_cast<int>(isqrt_u64(n)))) {
        alpha_cache_.reserve(1 << 20U);
    }

    u128 solve() { return dfs(1ULL, -1); }

  private:
    u64 n_;
    std::vector<int> primes_;
    std::unordered_map<u64, u64> alpha_cache_;

    u64 alpha(u64 n) {
        if (n <= 1ULL) return n;
        auto it = alpha_cache_.find(n);
        if (it != alpha_cache_.end()) return it->second;

        const u64 r = isqrt_u64(n);
        u128 s = 0;

        for (u64 i = 1; i <= r;) {
            const u64 q = n / i;
            const u64 j = std::min<u64>(r, n / q);
            s += static_cast<u128>(q) * static_cast<u128>(j - i + 1ULL);
            i = j + 1ULL;
        }

        const u64 ans = static_cast<u64>(2 * s - static_cast<u128>(r) * r);
        alpha_cache_.emplace(n, ans);
        return ans;
    }

    u128 dfs(u64 a, int i) {
        if (a > n_) return 0;

        u128 res = static_cast<u128>(alpha(n_ / a));

        for (int j = i + 1; j < static_cast<int>(primes_.size()); ++j) {
            const u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(j)]);
            if (a > n_ / p / p) break;

            u64 pp = p * p;
            u128 w = 1;

            while (a <= n_ / pp) {
                res += w * dfs(a * pp, j);
                if (pp > n_ / p) break;
                pp *= p;
                w <<= 1U;
            }
        }
        return res;
    }
};

u64 brute_sum(int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) continue;
        spf[static_cast<std::size_t>(i)] = i;
        if (1LL * i * i <= n) {
            for (int j = i * i; j <= n; j += i) {
                if (spf[static_cast<std::size_t>(j)] == 0) spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }

    u64 sum = 1;
    for (int x = 2; x <= n; ++x) {
        int t = x;
        u64 term = 1;
        while (t > 1) {
            int p = spf[static_cast<std::size_t>(t)];
            while (t % p == 0) {
                t /= p;
                term <<= 1ULL;
            }
        }
        sum += term;
    }
    return sum;
}

void run_validations() {
    assert(static_cast<u64>(Solver708(100ULL).solve()) == brute_sum(100));
    assert(static_cast<u64>(Solver708(50'000ULL).solve()) == brute_sum(50'000));
    assert(static_cast<u64>(Solver708(100'000'000ULL).solve()) == 9'613'563'919ULL);
}

}  // namespace

int main(int argc, char** argv) {
    bool validate = true;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--no-validate") validate = false;
    }

    if (validate) run_validations();

    constexpr u64 N = 100'000'000'000'000ULL;
    std::cout << to_string_u128(Solver708(N).solve()) << '\n';
    return 0;
}
