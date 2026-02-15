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
using i64 = std::int64_t;

constexpr u64 kDefaultN = 100'000'000'000ULL;
constexpr int kDefaultM = 100'000'000;

struct Options {
    u64 n = kDefaultN;
    int m = kDefaultM;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

std::string to_string_u128(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(v % 10)));
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }
    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }
    value = parsed;
    return true;
}

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<int>::max())) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--m=", options.m)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

std::vector<int> sieve_primes(const int n) {
    std::vector<bool> sieve(static_cast<std::size_t>(n) + 1, true);
    if (n >= 0) sieve[0] = false;
    if (n >= 1) sieve[1] = false;
    for (int p = 2; static_cast<i64>(p) * p <= n; ++p) {
        if (!sieve[static_cast<std::size_t>(p)]) continue;
        for (int q = p * p; q <= n; q += p) {
            sieve[static_cast<std::size_t>(q)] = false;
        }
    }
    std::vector<int> primes;
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (sieve[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    u64 x = a % mod;
    while (e > 0) {
        if (e & 1ULL) r = (r * x) % mod;
        x = (x * x) % mod;
        e >>= 1ULL;
    }
    return r;
}

u64 primitive_root_of_unity(u64 p, u64 d) {
    if (d == 1) return 1;
    for (u64 a = 2; ; ++a) {
        u64 h = mod_pow(a, (p - 1) / d, p);
        u64 z = 1;
        bool ok = true;
        for (u64 i = 1; i < d; ++i) {
            z = (z * h) % p;
            if (z == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return h;
    }
}

struct SumTask {
    const std::vector<int>* primes = nullptr;
    std::size_t start = 0;
    std::size_t end = 0;
    u64 n = 0;
    u128 sum = 0;
};

static void* sum_worker(void* arg) {
    auto* task = static_cast<SumTask*>(arg);
    const auto& primes = *task->primes;
    const u64 n = task->n;
    u128 local = 0;
    for (std::size_t i = task->start; i < task->end; ++i) {
        const u64 p = static_cast<u64>(primes[i]);
        const u64 d = gcd_u64(15, p - 1);
        const u64 r = primitive_root_of_unity(p, d);
        u64 h = p - 1;
        for (u64 j = 0; j < d; ++j) {
            h = (h * r) % p;
            local += static_cast<u128>((n + p - h) / p) * static_cast<u128>(p);
        }
    }
    task->sum = local;
    return nullptr;
}

unsigned thread_count(bool allow_multithreading, unsigned requested_threads) {
    if (!allow_multithreading) {
        return 1U;
    }
    unsigned threads = requested_threads;
    if (threads == 0U) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        if (n < 1) n = 1;
        threads = static_cast<unsigned>(n);
    }
    return std::max(1U, threads);
}

u128 sum_for(u64 n, int m, bool allow_multithreading, unsigned requested_threads) {
    const auto primes = sieve_primes(m);
    if (primes.empty()) {
        return 0;
    }

    const unsigned threads = thread_count(allow_multithreading, requested_threads);
    if (threads < 2U || primes.size() < 20000U) {
        u128 total = 0;
        for (int p_int : primes) {
            const u64 p = static_cast<u64>(p_int);
            const u64 d = gcd_u64(15, p - 1);
            const u64 r = primitive_root_of_unity(p, d);
            u64 h = p - 1;
            for (u64 j = 0; j < d; ++j) {
                h = (h * r) % p;
                total += static_cast<u128>((n + p - h) / p) * static_cast<u128>(p);
            }
        }
        return total;
    }

    const std::size_t count = primes.size();
    const unsigned active = std::min<unsigned>(threads, static_cast<unsigned>(count));
    const std::size_t block = (count + active - 1) / active;
    std::vector<SumTask> tasks(active);
    std::vector<pthread_t> workers(active);

    for (unsigned t = 0; t < active; ++t) {
        const std::size_t start = static_cast<std::size_t>(t) * block;
        const std::size_t end = std::min(count, start + block);
        tasks[t] = {&primes, start, end, n, 0};
        pthread_create(&workers[t], nullptr, sum_worker, &tasks[t]);
    }

    u128 total = 0;
    for (unsigned t = 0; t < active; ++t) {
        pthread_join(workers[t], nullptr);
        total += tasks[t].sum;
    }
    return total;
}

u64 s_single_bruteforce(u64 n, int m) {
    const auto primes = sieve_primes(m);
    u64 sum = 0;
    for (int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        const u64 r = mod_pow(n % p, 15ULL, p);
        if ((r + 1ULL) % p == 0ULL) sum += p;
    }
    return sum;
}

u128 brute_small(u64 n_limit, int m) {
    u128 s = 0;
    for (u64 n = 1; n <= n_limit; ++n) {
        s += s_single_bruteforce(n, m);
    }
    return s;
}

bool run_checkpoints() {
    if (s_single_bruteforce(2ULL, 10) != 3ULL) return false;
    if (s_single_bruteforce(2ULL, 1000) != 345ULL) return false;
    if (s_single_bruteforce(10ULL, 100) != 31ULL) return false;
    if (s_single_bruteforce(10ULL, 1000) != 483ULL) return false;
    if (sum_for(200ULL, 1000, false, 1) != brute_small(200ULL, 1000)) return false;
    if (sum_for(1000ULL, 5000, false, 1) != brute_small(1000ULL, 5000)) return false;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        std::cerr << "Checkpoint failed\n";
        return 2;
    }

    const u128 ans = sum_for(options.n, options.m, options.allow_multithreading, options.requested_threads);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
