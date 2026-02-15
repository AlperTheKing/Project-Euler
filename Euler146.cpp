#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 150000000ULL;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        u64 thread_count_value = 0;
        if (parse_u64_after_prefix(arg, "--threads=", thread_count_value)) {
            if (thread_count_value > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
                std::cerr << "Thread count too large: " << thread_count_value << '\n';
                return false;
            }
            options.requested_threads = static_cast<unsigned>(thread_count_value);
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 10;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    base %= mod;
    u64 result = 1 % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

bool is_prime(const u64 n) {
    if (n < 2) {
        return false;
    }
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0ULL) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    static constexpr std::array<u64, 7> kBases{
        2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};

    for (u64 a : kBases) {
        if (a % n == 0ULL) {
            continue;
        }
        u64 x = pow_mod(a, d, n);
        if (x == 1ULL || x == n - 1) {
            continue;
        }
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                witness = false;
                break;
            }
        }
        if (witness) {
            return false;
        }
    }
    return true;
}

std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; static_cast<long long>(p) * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) continue;
        for (int q = p * p; q <= n; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

struct ModFilter {
    int p = 0;
    int step = 0;
    std::vector<std::uint8_t> bad;
};

std::vector<ModFilter> build_filters(int limit) {
    static constexpr std::array<u64, 6> kGood{1ULL, 3ULL, 7ULL, 9ULL, 13ULL, 27ULL};
    const auto primes = sieve_primes(limit);
    std::vector<ModFilter> filters;
    filters.reserve(primes.size());
    for (int p : primes) {
        if (p == 2 || p == 5) {
            continue;
        }
        ModFilter f;
        f.p = p;
        f.step = 10 % p;
        f.bad.assign(static_cast<std::size_t>(p), 0);
        for (int r = 0; r < p; ++r) {
            const int r2 = static_cast<int>((1LL * r * r) % p);
            bool bad = false;
            for (u64 k : kGood) {
                if ((r2 + static_cast<int>(k % p)) % p == 0) {
                    bad = true;
                    break;
                }
            }
            f.bad[static_cast<std::size_t>(r)] = bad ? 1 : 0;
        }
        filters.push_back(std::move(f));
    }
    return filters;
}

unsigned thread_count(bool allow_multithreading, unsigned requested) {
    if (!allow_multithreading) {
        return 1U;
    }
    if (requested != 0U) {
        return requested;
    }
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    return static_cast<unsigned>(n);
}

struct Task {
    u64 start = 0;
    u64 end = 0;
    const std::vector<ModFilter>* filters = nullptr;
    u64 sum = 0;
};

u64 first_with_mod(u64 start, const u64 mod, const u64 residue) {
    const u64 current = start % mod;
    if (current <= residue) {
        return start + (residue - current);
    }
    return start + (mod - (current - residue));
}

void run_sequence(u64 start, const u64 end, const u64 residue, const std::vector<ModFilter>& filters, u64& sum) {
    static constexpr std::array<u64, 6> kGood{1ULL, 3ULL, 7ULL, 9ULL, 13ULL, 27ULL};
    static constexpr std::array<u64, 8> kBad{5ULL, 11ULL, 15ULL, 17ULL, 19ULL, 21ULL, 23ULL, 25ULL};

    constexpr u64 filter_threshold = 1000ULL;
    constexpr u64 step = 70ULL;

    u64 n = first_with_mod(start, step, residue);
    if (n >= end) {
        return;
    }

    u64 square = n * n;

    std::vector<int> residues(filters.size(), 0);
    std::vector<int> step70(filters.size(), 0);
    for (std::size_t i = 0; i < filters.size(); ++i) {
        const int p = filters[i].p;
        residues[i] = static_cast<int>(n % static_cast<u64>(p));
        step70[i] = (filters[i].step * 7) % p;
    }

    while (n < end) {
        if (square % 3ULL != 0ULL && square % 13ULL != 0ULL) {
            bool ok = true;
            if (n >= filter_threshold) {
                for (std::size_t i = 0; i < filters.size(); ++i) {
                    if (filters[i].bad[static_cast<std::size_t>(residues[i])]) {
                        ok = false;
                        break;
                    }
                }
            }
            if (ok) {
                for (u64 k : kGood) {
                    if (!is_prime(square + k)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    for (u64 k : kBad) {
                        if (is_prime(square + k)) {
                            ok = false;
                            break;
                        }
                    }
                }
                if (ok) {
                    sum += n;
                }
            }
        }

        n += step;
        square += 140ULL * (n - step) + 4900ULL;
        for (std::size_t i = 0; i < filters.size(); ++i) {
            const int p = filters[i].p;
            int next = residues[i] + step70[i];
            if (next >= p) next -= p;
            residues[i] = next;
        }
    }
}

static void* worker_fn(void* arg) {
    auto* task = static_cast<Task*>(arg);
    const auto& filters = *task->filters;

    u64 sum = 0;
    run_sequence(task->start, task->end, 10ULL, filters, sum);
    run_sequence(task->start, task->end, 60ULL, filters, sum);
    task->sum = sum;
    return nullptr;
}

u64 solve(const u64 limit, const std::vector<ModFilter>& filters, bool allow_multithreading, unsigned requested_threads) {
    const unsigned threads = thread_count(allow_multithreading, requested_threads);
    if (threads <= 1) {
        Task single{10, limit, &filters, 0};
        worker_fn(&single);
        return single.sum;
    }

    const u64 start = 10;
    const u64 count = (limit - start + 9) / 10;
    const u64 block = (count + threads - 1) / threads;

    std::vector<Task> tasks(threads);
    std::vector<pthread_t> workers(threads);
    for (unsigned t = 0; t < threads; ++t) {
        const u64 idx_start = static_cast<u64>(t) * block;
        const u64 idx_end = std::min(count, idx_start + block);
        const u64 n_start = start + idx_start * 10ULL;
        const u64 n_end = start + idx_end * 10ULL;
        tasks[t] = {n_start, n_end, &filters, 0};
        pthread_create(&workers[t], nullptr, worker_fn, &tasks[t]);
    }
    u64 sum = 0;
    for (unsigned t = 0; t < threads; ++t) {
        pthread_join(workers[t], nullptr);
        sum += tasks[t].sum;
    }
    return sum;
}

bool run_checkpoints() {
    if (!is_prime(2ULL) || !is_prime(3ULL) || is_prime(1ULL) || is_prime(21ULL)) {
        std::cerr << "Checkpoint failed for primality tester" << '\n';
        return false;
    }
    const auto filters = build_filters(2000);
    if (solve(1000000ULL, filters, false, 1) != 1242490ULL) {
        std::cerr << "Checkpoint failed for limit 1,000,000" << '\n';
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
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    const auto filters = build_filters(2000);
    std::cout << solve(options.limit, filters, options.allow_multithreading, options.requested_threads) << '\n';
    return 0;
}
