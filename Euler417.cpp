#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;
using i64 = std::int64_t;

struct Options {
    int limit = 100000000;
    int threads = 0;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 3 && options.threads >= 0;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0) {
        if (e & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        e >>= 1ULL;
    }
    return result;
}

int choose_thread_count(int requested, std::size_t work_items) {
    if (work_items <= 1) {
        return 1;
    }
    int threads = requested;
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
    }
    if (threads <= 0) {
        threads = 4;
    }
    if (threads > static_cast<int>(work_items)) {
        threads = static_cast<int>(work_items);
    }
    if (threads < 1) {
        threads = 1;
    }
    return threads;
}

std::vector<std::uint16_t> build_odd_spf(const int limit) {
    std::vector<std::uint16_t> spf(static_cast<std::size_t>(limit / 2 + 1), 0);
    const int root = static_cast<int>(std::sqrt(static_cast<double>(limit)));

    for (int i = 3; i <= root; i += 2) {
        if (spf[static_cast<std::size_t>(i >> 1)] != 0) {
            continue;
        }
        const i64 step = 2LL * i;
        for (i64 j = 1LL * i * i; j <= limit; j += step) {
            std::uint16_t& cell = spf[static_cast<std::size_t>(j >> 1)];
            if (cell == 0) {
                cell = static_cast<std::uint16_t>(i);
            }
        }
    }
    return spf;
}

std::vector<int> collect_primes(const int limit, const std::vector<std::uint16_t>& spf) {
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / 12));
    for (int p = 3; p <= limit; p += 2) {
        if (p != 5 && spf[static_cast<std::size_t>(p >> 1)] == 0) {
            primes.push_back(p);
        }
    }
    return primes;
}

u32 multiplicative_order_prime(const int p, const std::vector<std::uint16_t>& spf) {
    u32 ord = static_cast<u32>(p - 1);
    u32 rem = ord;

    if ((rem & 1U) == 0U) {
        while ((rem & 1U) == 0U) {
            rem >>= 1U;
        }
        while ((ord & 1U) == 0U) {
            const u32 cand = ord >> 1U;
            if (mod_pow(10, cand, static_cast<u64>(p)) == 1ULL) {
                ord = cand;
            } else {
                break;
            }
        }
    }

    while (rem > 1U) {
        u32 q = static_cast<u32>(spf[static_cast<std::size_t>(rem >> 1U)]);
        if (q == 0U) {
            q = rem;
        }
        while (rem % q == 0U) {
            rem /= q;
        }
        while (ord % q == 0U) {
            const u32 cand = ord / q;
            if (mod_pow(10, cand, static_cast<u64>(p)) == 1ULL) {
                ord = cand;
            } else {
                break;
            }
        }
    }
    return ord;
}

std::vector<u32> compute_prime_orders_parallel(const std::vector<int>& primes,
                                               const std::vector<std::uint16_t>& spf,
                                               const int limit,
                                               const int requested_threads) {
    std::vector<u32> ord_prime(static_cast<std::size_t>(limit / 2 + 1), 0U);
    ord_prime[0] = 1U;  // order(10 mod 1)

    const int thread_count = choose_thread_count(requested_threads, primes.size());
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        const std::size_t begin = primes.size() * static_cast<std::size_t>(t) /
                                  static_cast<std::size_t>(thread_count);
        const std::size_t end = primes.size() * static_cast<std::size_t>(t + 1) /
                                static_cast<std::size_t>(thread_count);
        workers.emplace_back([&, begin, end]() {
            for (std::size_t i = begin; i < end; ++i) {
                const int p = primes[i];
                ord_prime[static_cast<std::size_t>(p >> 1)] = multiplicative_order_prime(p, spf);
            }
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }
    return ord_prime;
}

std::vector<std::vector<u32>> build_small_prime_power_orders(
    const int limit, const std::vector<int>& primes, const std::vector<u32>& ord_prime) {
    const int root = static_cast<int>(std::sqrt(static_cast<double>(limit)));
    std::vector<std::vector<u32>> prime_power_orders(static_cast<std::size_t>(root + 1));

    for (const int p : primes) {
        if (p > root) {
            break;
        }
        auto& table = prime_power_orders[static_cast<std::size_t>(p)];
        int emax = 1;
        u64 pp = static_cast<u64>(p);
        while (pp <= static_cast<u64>(limit / p)) {
            pp *= static_cast<u64>(p);
            ++emax;
        }

        table.assign(static_cast<std::size_t>(emax + 1), 0U);
        table[1] = ord_prime[static_cast<std::size_t>(p >> 1)];

        pp = static_cast<u64>(p);
        for (int e = 2; e <= emax; ++e) {
            pp *= static_cast<u64>(p);
            const u32 prev = table[static_cast<std::size_t>(e - 1)];
            u32 cur = prev;
            if (mod_pow(10, prev, pp) != 1ULL) {
                cur = static_cast<u32>(prev * static_cast<u32>(p));
            }
            table[static_cast<std::size_t>(e)] = cur;
        }
    }
    return prime_power_orders;
}

std::vector<int> build_smooth_numbers(const int limit) {
    std::vector<int> smooth;
    for (u64 p2 = 1ULL; p2 <= static_cast<u64>(limit); p2 *= 2ULL) {
        for (u64 p5 = 1ULL; p2 * p5 <= static_cast<u64>(limit); p5 *= 5ULL) {
            smooth.push_back(static_cast<int>(p2 * p5));
        }
    }
    std::sort(smooth.begin(), smooth.end());
    smooth.erase(std::unique(smooth.begin(), smooth.end()), smooth.end());
    return smooth;
}

u32 multiplicative_order_from_factorization(
    int n, const std::vector<std::uint16_t>& spf, const std::vector<u32>& ord_prime,
    const std::vector<std::vector<u32>>& prime_power_orders) {
    u32 ord = 1U;
    int x = n;
    while (x > 1) {
        int p = static_cast<int>(spf[static_cast<std::size_t>(x >> 1)]);
        if (p == 0) {
            p = x;
        }
        int e = 1;
        int pp = p;
        x /= p;
        while (x % p == 0) {
            x /= p;
            pp *= p;
            ++e;
        }

        u32 part = 0U;
        if (e == 1) {
            part = ord_prime[static_cast<std::size_t>(p >> 1)];
        } else if (p < static_cast<int>(prime_power_orders.size()) &&
                   e < static_cast<int>(prime_power_orders[static_cast<std::size_t>(p)].size())) {
            part = prime_power_orders[static_cast<std::size_t>(p)][static_cast<std::size_t>(e)];
        } else {
            part = ord_prime[static_cast<std::size_t>(p >> 1)];
            u64 mod = static_cast<u64>(p);
            for (int k = 2; k <= e; ++k) {
                mod *= static_cast<u64>(p);
                if (mod_pow(10, part, mod) != 1ULL) {
                    part = static_cast<u32>(part * static_cast<u32>(p));
                }
            }
        }

        const u32 g = std::gcd(ord, part);
        ord = static_cast<u32>((static_cast<u64>(ord) / g) * part);
    }
    return ord;
}

u128 solve_fast(const int limit, const int requested_threads) {
    const std::vector<std::uint16_t> spf = build_odd_spf(limit);
    const std::vector<int> primes = collect_primes(limit, spf);
    const std::vector<u32> ord_prime =
        compute_prime_orders_parallel(primes, spf, limit, requested_threads);
    const std::vector<std::vector<u32>> prime_power_orders =
        build_small_prime_power_orders(limit, primes, ord_prime);
    const std::vector<int> smooth = build_smooth_numbers(limit);

    const std::size_t odd_count =
        static_cast<std::size_t>((static_cast<u64>(limit) - 3ULL) / 2ULL + 1ULL);
    const int thread_count = choose_thread_count(requested_threads, odd_count);
    std::vector<u128> partial(static_cast<std::size_t>(thread_count), 0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        const std::size_t begin = odd_count * static_cast<std::size_t>(t) /
                                  static_cast<std::size_t>(thread_count);
        const std::size_t end = odd_count * static_cast<std::size_t>(t + 1) /
                                static_cast<std::size_t>(thread_count);
        workers.emplace_back([&, begin, end, t]() {
            u128 local = 0;
            for (std::size_t idx = begin; idx < end; ++idx) {
                const int n = static_cast<int>(3 + 2 * idx);
                if (n % 5 == 0) {
                    continue;
                }
                const u32 ord =
                    multiplicative_order_from_factorization(n, spf, ord_prime, prime_power_orders);
                const int quotient = limit / n;
                const int multiplicity =
                    static_cast<int>(std::upper_bound(smooth.begin(), smooth.end(), quotient) -
                                     smooth.begin());
                local += static_cast<u128>(ord) * static_cast<u128>(multiplicity);
            }
            partial[static_cast<std::size_t>(t)] = local;
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }

    u128 total = 0;
    for (const u128 value : partial) {
        total += value;
    }
    return total;
}

u64 recurring_cycle_length_bruteforce(int n) {
    while ((n % 2) == 0) {
        n /= 2;
    }
    while ((n % 5) == 0) {
        n /= 5;
    }
    if (n == 1) {
        return 0;
    }

    int r = 10 % n;
    u64 len = 1;
    while (r != 1) {
        r = static_cast<int>((static_cast<i64>(r) * 10LL) % n);
        ++len;
    }
    return len;
}

u128 solve_bruteforce(const int limit) {
    u128 total = 0;
    for (int n = 3; n <= limit; ++n) {
        total += recurring_cycle_length_bruteforce(n);
    }
    return total;
}

bool run_checkpoints(const int requested_threads) {
    if (solve_bruteforce(10) != 9) {
        std::cerr << "Checkpoint failed: sum_{n=3..10} L(n)\n";
        return false;
    }

    if (solve_fast(1000, requested_threads) != solve_bruteforce(1000)) {
        std::cerr << "Checkpoint failed: fast solver vs brute force at limit=1000\n";
        return false;
    }

    if (solve_fast(1000000, requested_threads) != static_cast<u128>(55535191115ULL)) {
        std::cerr << "Checkpoint failed: sum_{n=3..1000000} L(n)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints(options.threads)) {
        return 2;
    }

    const u128 answer = solve_fast(options.limit, options.threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
