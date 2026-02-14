#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct GaussInt {
    i128 x;
    i128 y;
};

static inline GaussInt multiply(const GaussInt& a, const GaussInt& b) {
    return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
}

static inline GaussInt conjugate(const GaussInt& a) {
    return {a.x, -a.y};
}

static u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) * (r + 1ULL) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

static std::vector<int> prime_list(const int limit) {
    std::vector<char> composite(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (composite[static_cast<std::size_t>(i)]) continue;
        primes.push_back(i);
        if (1LL * i * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[static_cast<std::size_t>(j)] = 1;
            }
        }
    }
    return primes;
}

static std::pair<int, int> sum_of_two_squares_prime(const int p) {
    for (int x = 1; 1LL * x * x <= p; ++x) {
        const int y2 = p - x * x;
        const u64 y = isqrt_u64(static_cast<u64>(y2));
        if (static_cast<int>(y * y) == y2) return {x, static_cast<int>(y)};
    }
    return {0, 0};
}

static u64 factor_with_primes(u64 value,
                              const std::vector<int>& primes,
                              const int max_prime,
                              std::vector<std::pair<int, int>>& factors) {
    factors.clear();
    for (const int p : primes) {
        if (p > max_prime) break;
        if (value % static_cast<u64>(p) != 0ULL) continue;
        int e = 0;
        do {
            value /= static_cast<u64>(p);
            ++e;
        } while (value % static_cast<u64>(p) == 0ULL);
        factors.push_back({p, e});
    }
    return value;
}

static int count_from_factorization(const std::vector<std::pair<int, int>>& factors,
                                    const std::map<int, GaussInt>& gaussian_primes) {
    std::vector<GaussInt> reps;
    reps.push_back({1, 0});

    for (const auto& [prime, exp] : factors) {
        const GaussInt base = gaussian_primes.at(prime);
        std::vector<GaussInt> powers(static_cast<std::size_t>(exp + 1));
        powers[0] = {1, 0};
        for (int i = 1; i <= exp; ++i) {
            powers[static_cast<std::size_t>(i)] = multiply(powers[static_cast<std::size_t>(i - 1)], base);
        }

        std::vector<GaussInt> terms(static_cast<std::size_t>(exp + 1));
        for (int i = 0; i <= exp; ++i) {
            terms[static_cast<std::size_t>(i)] =
                multiply(powers[static_cast<std::size_t>(i)],
                         conjugate(powers[static_cast<std::size_t>(exp - i)]));
        }

        std::vector<GaussInt> next;
        next.reserve(reps.size() * static_cast<std::size_t>(exp + 1));
        for (const auto& r : reps) {
            for (const auto& t : terms) {
                next.push_back(multiply(r, t));
            }
        }
        reps.swap(next);
    }

    std::set<std::pair<u64, u64>> uniq;
    for (const auto& rep : reps) {
        i128 x = rep.x < 0 ? -rep.x : rep.x;
        i128 y = rep.y < 0 ? -rep.y : rep.y;
        if (x == 0 || y == 0) continue;

        u64 xx = static_cast<u64>(x);
        u64 yy = static_cast<u64>(y);
        if (xx < yy) std::swap(xx, yy);
        if (xx % 6ULL == 5ULL && yy % 6ULL == 5ULL) {
            uniq.insert({xx, yy});
        }
    }
    return static_cast<int>(uniq.size());
}

static int brute_count(const int n) {
    auto pent = [](const int x) -> u64 { return static_cast<u64>(x) * (3ULL * static_cast<u64>(x) - 1ULL) / 2ULL; };
    const u64 target = pent(n);
    int cnt = 0;
    for (int a = 1; a <= n; ++a) {
        for (int b = 1; b <= a; ++b) {
            if (pent(a) + pent(b) == target) ++cnt;
        }
    }
    return cnt;
}

static int count_exact_small_n(const int n, std::map<int, GaussInt>& gaussian_primes) {
    const u64 k = 6ULL * static_cast<u64>(n) - 1ULL;
    u64 value = k * k + 1ULL;

    std::vector<std::pair<int, int>> factors;
    for (u64 p = 2; p * p <= value; ++p) {
        if (value % p != 0ULL) continue;
        int e = 0;
        do {
            value /= p;
            ++e;
        } while (value % p == 0ULL);

        const int ip = static_cast<int>(p);
        factors.push_back({ip, e});
        if (ip == 2) continue;
        if (ip % 4 == 3) return 0;
        if (!gaussian_primes.count(ip)) {
            const auto [x, y] = sum_of_two_squares_prime(ip);
            gaussian_primes[ip] = {x, y};
        }
    }
    if (value > 1ULL) {
        const int p = static_cast<int>(value);
        factors.push_back({p, 1});
        if (p != 2 && p % 4 == 3) return 0;
        if (p == 2) {
            gaussian_primes[p] = {1, 1};
        } else if (!gaussian_primes.count(p)) {
            const auto [x, y] = sum_of_two_squares_prime(p);
            gaussian_primes[p] = {x, y};
        }
    }
    return count_from_factorization(factors, gaussian_primes);
}

int main() {
    constexpr int threshold = 100;
    constexpr int prime_limit = 400;

    const std::vector<int> primes = prime_list(prime_limit);
    std::map<int, GaussInt> gaussian_primes;
    gaussian_primes[2] = {1, 1};

    for (const int p : primes) {
        if (p == 2 || p % 4 == 3) continue;
        const auto [x, y] = sum_of_two_squares_prime(p);
        gaussian_primes[p] = {x, y};
    }

    assert(count_exact_small_n(49, gaussian_primes) == 2);
    assert(count_exact_small_n(268, gaussian_primes) == 3);
    for (int n = 1; n <= 120; ++n) {
        assert(count_exact_small_n(n, gaussian_primes) == brute_count(n));
    }

    std::vector<std::pair<int, int>> factors;
    factors.reserve(16);

    for (u64 n = 1;; ++n) {
        const u64 k = 6ULL * n - 1ULL;
        const u64 value = k * k + 1ULL;

        const u64 rem = factor_with_primes(value, primes, prime_limit, factors);
        if (rem != 1ULL) continue;

        int bound = 1;
        for (const auto& [_, e] : factors) {
            bound *= (e + 1);
            if (bound > 2 * (threshold + 1)) break;
        }
        if (bound / 2 <= threshold) continue;

        const int ways = count_from_factorization(factors, gaussian_primes);
        if (ways > threshold) {
            const u64 answer = n * (3ULL * n - 1ULL) / 2ULL;
            std::cout << answer << '\n';
            return 0;
        }
    }
}
