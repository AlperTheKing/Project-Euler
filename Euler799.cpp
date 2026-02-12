#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

static inline u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod(r, a, mod);
        }
        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    static constexpr u64 WITNESSES[] = {2ULL, 325ULL, 9'375ULL, 28'178ULL, 450'775ULL, 9'780'504ULL, 1'795'265'022ULL};
    for (u64 a : WITNESSES) {
        if (a % n == 0) {
            continue;
        }

        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) {
            continue;
        }

        bool composite = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) {
            return false;
        }
    }

    return true;
}

static u64 pollard_rho(u64 n, std::mt19937_64& rng) {
    if ((n & 1ULL) == 0ULL) {
        return 2ULL;
    }
    if (n % 3ULL == 0ULL) {
        return 3ULL;
    }

    std::uniform_int_distribution<u64> dist(2ULL, n - 2ULL);

    while (true) {
        const u64 c = dist(rng);
        u64 x = dist(rng);
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) {
            return (mul_mod(v, v, n) + c) % n;
        };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            const u64 diff = x > y ? x - y : y - x;
            d = std::gcd(diff, n);
        }

        if (d != n) {
            return d;
        }
    }
}

static void factor_rec(u64 n, std::map<u64, int>& out, std::mt19937_64& rng) {
    if (n == 1) {
        return;
    }
    if (is_prime(n)) {
        ++out[n];
        return;
    }

    const u64 d = pollard_rho(n, rng);
    factor_rec(d, out, rng);
    factor_rec(n / d, out, rng);
}

static u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while (static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > x) {
        --r;
    }
    return r;
}

static std::pair<i128, i128> gaussian_mul(const std::pair<i128, i128>& a, const std::pair<i128, i128>& b) {
    return {a.first * b.first - a.second * b.second, a.first * b.second + a.second * b.first};
}

static std::pair<u64, u64> prime_as_sum_of_two_squares(u64 p, std::unordered_map<u64, std::pair<u64, u64>>& cache) {
    const auto it = cache.find(p);
    if (it != cache.end()) {
        return it->second;
    }

    u64 a = 2;
    while (pow_mod(a, (p - 1ULL) / 2ULL, p) != p - 1ULL) {
        ++a;
    }
    const u64 r = pow_mod(a, (p - 1ULL) / 4ULL, p);

    i128 A = static_cast<i128>(p);
    i128 B = static_cast<i128>(r);
    while (B * B > static_cast<i128>(p)) {
        const i128 t = A % B;
        A = B;
        B = t;
    }

    const u64 x = static_cast<u64>(B);
    const u64 y2 = p - x * x;
    const u64 y = isqrt_u64(y2);
    assert(static_cast<u128>(y) * static_cast<u128>(y) == y2);

    cache[p] = {x, y};
    return {x, y};
}

static int count_representations_fast(u64 n,
                                      std::mt19937_64& rng,
                                      std::unordered_map<u64, std::pair<u64, u64>>& sumsq_cache,
                                      int threshold,
                                      int* out_divisor_product) {
    const u64 z = 6ULL * n - 1ULL;
    const u64 norm = z * z + 1ULL;

    std::map<u64, int> fac;
    factor_rec(norm, fac, rng);

    int divisor_product = 1;
    for (const auto& [p, e] : fac) {
        if (p == 2ULL) {
            continue;
        }
        divisor_product *= (e + 1);
    }
    *out_divisor_product = divisor_product;

    if (threshold >= 0 && divisor_product / 2 <= threshold) {
        return divisor_product / 2;
    }

    std::vector<std::pair<i128, i128>> states;
    states.push_back({1, 1});

    for (const auto& [p, e] : fac) {
        if (p == 2ULL) {
            continue;
        }

        const auto [a, b] = prime_as_sum_of_two_squares(p, sumsq_cache);
        const std::pair<i128, i128> pi = {static_cast<i128>(a), static_cast<i128>(b)};
        const std::pair<i128, i128> conj_pi = {static_cast<i128>(a), -static_cast<i128>(b)};

        std::vector<std::pair<i128, i128>> pow_pi(e + 1), pow_conj(e + 1);
        pow_pi[0] = {1, 0};
        pow_conj[0] = {1, 0};

        for (int i = 1; i <= e; ++i) {
            pow_pi[i] = gaussian_mul(pow_pi[i - 1], pi);
            pow_conj[i] = gaussian_mul(pow_conj[i - 1], conj_pi);
        }

        std::vector<std::pair<i128, i128>> next;
        next.reserve(states.size() * static_cast<std::size_t>(e + 1));

        for (const auto& g : states) {
            for (int k = 0; k <= e; ++k) {
                auto t = gaussian_mul(g, pow_pi[k]);
                t = gaussian_mul(t, pow_conj[e - k]);
                next.push_back(t);
            }
        }

        states.swap(next);
    }

    std::set<std::pair<u64, u64>> uniq;
    for (auto g : states) {
        i128 x = g.first;
        i128 y = g.second;
        if (x < 0) {
            x = -x;
        }
        if (y < 0) {
            y = -y;
        }
        if (x == 0 || y == 0) {
            continue;
        }

        u64 xx = static_cast<u64>(x);
        u64 yy = static_cast<u64>(y);
        if (xx < yy) {
            std::swap(xx, yy);
        }
        if (xx % 6ULL == 5ULL && yy % 6ULL == 5ULL) {
            uniq.insert({xx, yy});
        }
    }

    return static_cast<int>(uniq.size());
}

static int count_representations_bruteforce(int n) {
    std::vector<u64> pent(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        pent[i] = static_cast<u64>(i) * static_cast<u64>(3 * i - 1) / 2ULL;
    }

    int cnt = 0;
    for (int a = 1; a <= n; ++a) {
        for (int b = 1; b <= a; ++b) {
            const u64 s = pent[a] + pent[b];
            const u64 d = 24ULL * s + 1ULL;
            const u64 r = isqrt_u64(d);
            if (r * r != d) {
                continue;
            }
            if ((1ULL + r) % 6ULL != 0ULL) {
                continue;
            }
            if ((1ULL + r) / 6ULL == static_cast<u64>(n)) {
                ++cnt;
            }
        }
    }
    return cnt;
}

static u64 find_first_index_over_threshold(int threshold) {
    int thread_count = static_cast<int>(std::thread::hardware_concurrency());
    if (thread_count <= 0) {
        thread_count = 4;
    }

    const u64 block_size = 200'000ULL;
    u64 block_start = 1ULL;

    while (true) {
        const u64 block_end = block_start + block_size - 1ULL;
        std::atomic<u64> next(block_start);
        std::atomic<u64> best_in_block(std::numeric_limits<u64>::max());

        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(thread_count));

        for (int t = 0; t < thread_count; ++t) {
            workers.emplace_back([&, t]() {
                std::mt19937_64 rng(0x9E3779B97F4A7C15ULL ^ (static_cast<u64>(t) << 32) ^ block_start);
                std::unordered_map<u64, std::pair<u64, u64>> sumsq_cache;
                sumsq_cache.reserve(2048);

                while (true) {
                    const u64 n = next.fetch_add(1ULL, std::memory_order_relaxed);
                    if (n > block_end) {
                        break;
                    }

                    const u64 best_now = best_in_block.load(std::memory_order_relaxed);
                    if (n >= best_now) {
                        break;
                    }

                    int divisor_product = 0;
                    const int ways = count_representations_fast(n, rng, sumsq_cache, threshold, &divisor_product);
                    if (ways > threshold) {
                        u64 cur = best_in_block.load(std::memory_order_relaxed);
                        while (n < cur && !best_in_block.compare_exchange_weak(cur, n, std::memory_order_relaxed)) {
                        }
                    }
                }
            });
        }

        for (auto& th : workers) {
            th.join();
        }

        const u64 candidate = best_in_block.load(std::memory_order_relaxed);
        if (candidate != std::numeric_limits<u64>::max()) {
            return candidate;
        }

        block_start = block_end + 1ULL;
    }
}

int main() {
    std::mt19937_64 rng(0xC0FFEE123456789ULL);
    std::unordered_map<u64, std::pair<u64, u64>> sumsq_cache;
    sumsq_cache.reserve(1024);

    int div_prod = 0;
    assert(count_representations_fast(49ULL, rng, sumsq_cache, -1, &div_prod) == 2);
    assert(count_representations_fast(268ULL, rng, sumsq_cache, -1, &div_prod) == 3);

    for (int n = 1; n <= 120; ++n) {
        const int fast = count_representations_fast(static_cast<u64>(n), rng, sumsq_cache, -1, &div_prod);
        const int brute = count_representations_bruteforce(n);
        assert(fast == brute);
    }

    const u64 n = find_first_index_over_threshold(100);
    const u64 answer = n * (3ULL * n - 1ULL) / 2ULL;
    std::cout << answer << '\n';
    return 0;
}
