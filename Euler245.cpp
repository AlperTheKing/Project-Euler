#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = uint64_t;
using u128 = __uint128_t;

constexpr u64 kDefaultLimit = 200000000000ULL;

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

u64 icbrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) * (r + 1) <= x) ++r;
    while (r * r * r > x) --r;
    return r;
}

std::vector<int> sieve_primes(int max_n, std::vector<bool>* is_prime_out) {
    std::vector<bool> is_prime(max_n + 1, true);
    is_prime[0] = false;
    if (max_n >= 1) is_prime[1] = false;
    for (int i = 2; (int64_t)i * i <= max_n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= max_n; j += i) is_prime[j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= max_n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    if (is_prime_out) *is_prime_out = std::move(is_prime);
    return primes;
}

u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 result = 1 % mod;
    u64 base = a % mod;
    while (e > 0) {
        if (e & 1) result = mod_mul(result, base, mod);
        base = mod_mul(base, base, mod);
        e >>= 1;
    }
    return result;
}

bool is_prime_mr(u64 n) {
    if (n < 2) return false;
    static const u64 small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 p : small_primes) {
        if (n % p == 0) return n == p;
    }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }
    auto witness = [&](u64 a) {
        if (a % n == 0) return false;
        u64 x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) return false;
        for (int i = 1; i < s; ++i) {
            x = mod_mul(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };
    static const u64 bases[] = {2, 3, 5, 7, 11, 13};
    for (u64 a : bases) {
        if (witness(a)) return false;
    }
    return true;
}

bool is_prime(u64 n, const std::vector<bool>& is_prime_small) {
    if (n < is_prime_small.size()) return is_prime_small[n];
    return is_prime_mr(n);
}

u64 tonelli_shanks(u64 n, u64 p) {
    if (n == 0) return 0;
    if (p == 2) return 1;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);
    u64 q = p - 1;
    int s = 0;
    while ((q & 1) == 0) {
        q >>= 1;
        ++s;
    }
    u64 z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) ++z;
    u64 c = mod_pow(z, q, p);
    u64 x = mod_pow(n, (q + 1) / 2, p);
    u64 t = mod_pow(n, q, p);
    int m = s;
    while (t != 1) {
        int i = 1;
        u64 t2i = mod_mul(t, t, p);
        while (t2i != 1) {
            t2i = mod_mul(t2i, t2i, p);
            ++i;
        }
        u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        x = mod_mul(x, b, p);
        u64 b2 = mod_mul(b, b, p);
        t = mod_mul(t, b2, p);
        c = b2;
        m = i;
    }
    return x;
}

struct Factor {
    u64 prime = 0;
    int exp = 0;
};

void collect_semiprimes(u64 limit, int threads, const std::vector<int>& primes,
                        const std::vector<bool>& is_prime_small, std::vector<u64>& out) {
    u64 k_max = isqrt_u64(limit);
    std::vector<u64> rem(k_max + 1);
    std::vector<std::vector<Factor>> factors(k_max + 1);
    for (u64 k = 2; k <= k_max; ++k) {
        rem[k] = k * k - k + 1;
    }

    for (int p : primes) {
        if (p > static_cast<int>(k_max)) break;
        if (p == 2) continue;
        if (p == 3) {
            for (u64 k = 2; k <= k_max; k += 3) {
                if (rem[k] % 3 != 0) continue;
                int exp = 0;
                while (rem[k] % 3 == 0) {
                    rem[k] /= 3;
                    ++exp;
                }
                factors[k].push_back({3, exp});
            }
            continue;
        }
        if (p % 3 != 1) continue;

        u64 n = static_cast<u64>(p - 3);
        u64 root = tonelli_shanks(n, static_cast<u64>(p));
        u64 inv2 = (static_cast<u64>(p) + 1) / 2;
        u64 r1 = (1 + root) % p;
        r1 = (r1 * inv2) % p;
        u64 r2 = (1 + p - root) % p;
        r2 = (r2 * inv2) % p;

        auto sieve_root = [&](u64 r) {
            if (r < 2) {
                u64 add = (2 - r + p - 1) / p;
                r += add * static_cast<u64>(p);
            }
            for (u64 k = r; k <= k_max; k += p) {
                if (rem[k] % p != 0) continue;
                int exp = 0;
                while (rem[k] % p == 0) {
                    rem[k] /= p;
                    ++exp;
                }
                factors[k].push_back({static_cast<u64>(p), exp});
            }
        };

        sieve_root(r1);
        if (r2 != r1) sieve_root(r2);
    }

    for (u64 k = 2; k <= k_max; ++k) {
        if (rem[k] > 1) factors[k].push_back({rem[k], 1});
    }

    int use_threads = std::max(1, threads);
    std::vector<std::vector<u64>> locals(use_threads);
    auto worker = [&](int tid, u64 start, u64 end) {
        std::vector<u64> divisors;
        divisors.reserve(64);
        for (u64 k = start; k < end; ++k) {
            u64 K = k * k - k + 1;
            divisors.clear();
            divisors.push_back(1);
            for (const auto& f : factors[k]) {
                size_t size = divisors.size();
                u64 mult = 1;
                for (int e = 1; e <= f.exp; ++e) {
                    mult *= f.prime;
                    for (size_t i = 0; i < size; ++i) {
                        divisors.push_back(divisors[i] * mult);
                    }
                }
            }
            for (u64 a : divisors) {
                u64 b = K / a;
                if (a > b) continue;
                u64 p = k + a;
                u64 q = k + b;
                if (p == q) continue;
                if (static_cast<u128>(p) * q > limit) continue;
                if (is_prime(p, is_prime_small) && is_prime(q, is_prime_small)) {
                    locals[tid].push_back(p * q);
                }
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(use_threads);
    u64 total = k_max - 1;
    u64 chunk = (total + use_threads - 1) / use_threads;
    for (int t = 0; t < use_threads; ++t) {
        u64 start = 2 + t * chunk;
        u64 end = std::min<u64>(2 + (t + 1) * chunk, k_max + 1);
        if (start >= end) break;
        pool.emplace_back(worker, t, start, end);
    }
    for (auto& th : pool) th.join();

    for (auto& vec : locals) {
        out.insert(out.end(), vec.begin(), vec.end());
    }
}

void dfs_multi(int k, size_t idx, u64 m, u64 A, int last_prime, int depth, u64 limit,
               const std::vector<int>& primes, const std::vector<bool>& is_prime_small,
               std::vector<u64>& out) {
    u128 denom = static_cast<u128>(k) * A - static_cast<u128>(k - 1) * m;
    if (denom <= 0) return;
    u128 numer = static_cast<u128>(k) * A + 1;
    if (numer % denom == 0) {
        u64 p = static_cast<u64>(numer / denom);
        if (p > static_cast<u64>(last_prime) && static_cast<u128>(m) * p <= limit) {
            if (depth >= 2 && is_prime(p, is_prime_small)) {
                out.push_back(m * p);
            }
        }
    }

    for (size_t i = idx; i < primes.size(); ++i) {
        int q = primes[i];
        if (q <= k) continue;
        if (static_cast<u128>(m) * q * q > limit) break;
        dfs_multi(k, i + 1, m * static_cast<u64>(q), A * static_cast<u64>(q - 1), q,
                  depth + 1, limit, primes, is_prime_small, out);
    }
}

void collect_multi_primes(u64 limit, int threads, const std::vector<int>& primes,
                          const std::vector<bool>& is_prime_small, std::vector<u64>& out) {
    u64 k_root = icbrt_u64(limit);
    if (k_root == 0) return;
    int k_max = static_cast<int>(k_root) - 1;
    if (k_max < 2) return;

    int use_threads = std::max(1, threads);
    std::vector<std::vector<u64>> locals(use_threads);
    auto worker = [&](int tid, int start_k, int end_k) {
        for (int k = start_k; k < end_k; ++k) {
            auto it = std::upper_bound(primes.begin(), primes.end(), k);
            size_t idx = static_cast<size_t>(it - primes.begin());
            dfs_multi(k, idx, 1, 1, 1, 0, limit, primes, is_prime_small, locals[tid]);
        }
    };

    int total = k_max - 1;
    int chunk = (total + use_threads - 1) / use_threads;
    std::vector<std::thread> pool;
    for (int t = 0; t < use_threads; ++t) {
        int start_k = 2 + t * chunk;
        int end_k = std::min(2 + (t + 1) * chunk, k_max + 1);
        if (start_k >= end_k) break;
        pool.emplace_back(worker, t, start_k, end_k);
    }
    for (auto& th : pool) th.join();

    for (auto& vec : locals) {
        out.insert(out.end(), vec.begin(), vec.end());
    }
}

u64 solve_fast(u64 limit, int threads) {
    u64 k_max = isqrt_u64(limit);
    std::vector<bool> is_prime_small;
    std::vector<int> primes = sieve_primes(static_cast<int>(k_max), &is_prime_small);
    std::vector<u64> results;
    results.reserve(1024);
    collect_semiprimes(limit, threads, primes, is_prime_small, results);
    collect_multi_primes(limit, threads, primes, is_prime_small, results);
    std::sort(results.begin(), results.end());
    results.erase(std::unique(results.begin(), results.end()), results.end());
    u128 sum = 0;
    for (u64 n : results) sum += n;
    return static_cast<u64>(sum);
}

u64 brute_sum(u64 limit) {
    std::vector<u64> phi(limit + 1);
    for (u64 i = 0; i <= limit; ++i) phi[i] = i;
    for (u64 i = 2; i <= limit; ++i) {
        if (phi[i] != i) continue;
        for (u64 j = i; j <= limit; j += i) {
            phi[j] -= phi[j] / i;
        }
    }
    u64 sum = 0;
    for (u64 n = 4; n <= limit; ++n) {
        u64 c = n - phi[n];
        if ((n - 1) % c != 0) continue;
        if (phi[n] == n - 1) continue;
        sum += n;
    }
    return sum;
}

void run_validation() {
    struct Check {
        u64 limit;
        u64 expected;
    };
    const Check checks[] = {
        {1000, 1594},
        {100000, 595394},
        {1000000, 11149065},
    };
    for (const auto& chk : checks) {
        u64 brute = brute_sum(chk.limit);
        if (brute != chk.expected) {
            std::cerr << "Validation failed (brute) for N=" << chk.limit
                      << ": got " << brute << " expected " << chk.expected << "\n";
            std::exit(1);
        }
        u64 fast = solve_fast(chk.limit, 1);
        if (fast != chk.expected) {
            std::cerr << "Validation failed (fast) for N=" << chk.limit
                      << ": got " << fast << " expected " << chk.expected << "\n";
            std::exit(1);
        }
    }
    std::cerr << "Validation checkpoints passed.\n";
}

void usage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " [-n LIMIT] [-t THREADS] [--no-validate]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    u64 limit = kDefaultLimit;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-n" && i + 1 < argc) {
            limit = std::stoull(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--validate") {
            validate = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (validate) run_validation();

    u64 answer = solve_fast(limit, threads);
    std::cout << answer << "\n";
    return 0;
}
