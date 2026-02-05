#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

namespace {

constexpr u64 kDefaultLimit = 1000000000000000000ULL;  // 1e18

u64 gcd_u64(u64 a, u64 b) {
    while (b) {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u128 gcd_u128(u128 a, u128 b) {
    while (b) {
        u128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u64 mul_mod_u64(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 pow_mod_u64(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) {
        if (d & 1) r = mul_mod_u64(r, a, mod);
        a = mul_mod_u64(a, a, mod);
        d >>= 1;
    }
    return r;
}

bool is_prime_u64(u64 n) {
    if (n < 2) return false;
    static u64 small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    u64 d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](u64 a) -> bool {
        if (a % n == 0) return false;
        u64 x = pow_mod_u64(a, d, n);
        if (x == 1 || x == n - 1) return false;
        for (u64 i = 1; i < s; ++i) {
            x = mul_mod_u64(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };

    static u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (witness(a)) return false;
    }
    return true;
}

u64 splitmix64(u64& x) {
    u64 z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

thread_local u64 rng_state =
    static_cast<u64>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) ^
    static_cast<u64>(std::hash<std::thread::id>()(std::this_thread::get_id()));

u64 rand_u64(u64 lo, u64 hi) {
    u64 r = splitmix64(rng_state);
    return lo + (hi > lo ? (r % (hi - lo + 1)) : 0);
}

u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    if (n % 3ULL == 0) return 3;

    while (true) {
        u64 c = rand_u64(1, n - 1);
        u64 x = rand_u64(0, n - 1);
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) { return (mul_mod_u64(v, v, n) + c) % n; };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            u64 diff = (x > y) ? (x - y) : (y - x);
            d = gcd_u64(diff, n);
        }
        if (d != n) return d;
    }
}

u64 min_factor_u64(u64 n);

u64 min_factor_rec(u64 n) {
    if (n == 1) return 1;
    if (is_prime_u64(n)) return n;
    u64 d = pollard_rho(n);
    u64 a = min_factor_rec(d);
    u64 b = min_factor_rec(n / d);
    return std::min(a, b);
}

u64 min_factor_u64(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    return min_factor_rec(n);
}

struct Solver {
    u64 limit;
    u64 target_num;  // A (odd), target ratio is A/2.
    u128 sum = 0;
    std::unordered_map<u64, u64> min_factor_cache;

    explicit Solver(u64 limit_, u64 target_num_) : limit(limit_), target_num(target_num_) {
        min_factor_cache.reserve(1 << 16);
    }

    u64 get_min_factor(u64 n) {
        auto it = min_factor_cache.find(n);
        if (it != min_factor_cache.end()) return it->second;
        u64 res = min_factor_u64(n);
        min_factor_cache.emplace(n, res);
        return res;
    }

    void dfs(u64 n, u64 rn, u64 rs) {
        if (static_cast<u128>(n) * rs > limit) return;
        if (rn == rs) {
            sum += n;
            return;
        }
        if (rn < rs) return;  // overshot target
        if (rs == 1) return;  // no forced prime left; with n even this can't finish within limit

        u64 p = get_min_factor(rs);
        if (n % p == 0) return;  // p already used; would have been chosen with larger exponent earlier

        u64 tmp = rs;
        int e = 0;
        while (tmp % p == 0) {
            tmp /= p;
            ++e;
        }

        u128 p_pow = 1;
        u128 sigma = 1;
        for (int i = 1; ; ++i) {
            if (p_pow > limit / p) break;
            p_pow *= p;
            sigma += p_pow;

            if (i < e) continue;
            u128 n2 = static_cast<u128>(n) * p_pow;
            if (n2 > limit) break;

            u128 rn2 = static_cast<u128>(rn) * p_pow;
            u128 rs2 = static_cast<u128>(rs) * sigma;
            u128 g = gcd_u128(rn2, rs2);
            rn2 /= g;
            rs2 /= g;

            if (rn2 < rs2) break;

            if (n2 * rs2 > limit) continue;
            dfs(static_cast<u64>(n2), static_cast<u64>(rn2), static_cast<u64>(rs2));
        }
    }

    u128 run() {
        dfs(1ULL, target_num, 2ULL);
        return sum;
    }
};

std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 solve(u64 limit, int threads) {
    std::vector<u64> targets;
    for (u64 a = 3; a <= 13; a += 2) targets.push_back(a);

    if (threads <= 1 || targets.size() == 1) {
        u128 total = 0;
        for (u64 a : targets) {
            Solver solver(limit, a);
            total += solver.run();
        }
        return total;
    }

    int worker_count = std::min<int>(threads, static_cast<int>(targets.size()));
    std::vector<u128> partial(targets.size(), 0);
    std::atomic<std::size_t> idx{0};

    auto worker = [&]() {
        while (true) {
            std::size_t i = idx.fetch_add(1);
            if (i >= targets.size()) return;
            Solver solver(limit, targets[i]);
            partial[i] = solver.run();
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(worker_count);
    for (int i = 0; i < worker_count; ++i) {
        pool.emplace_back(worker);
    }
    for (auto& t : pool) t.join();

    u128 total = 0;
    for (u128 v : partial) total += v;
    return total;
}

void usage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " [-n LIMIT] [-t THREADS]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    u64 limit = kDefaultLimit;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-n" && i + 1 < argc) {
            limit = std::stoull(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            threads = std::max(1, std::stoi(argv[++i]));
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    u128 ans = solve(limit, threads);
    std::cout << to_string_u128(ans) << "\n";
    return 0;
}
