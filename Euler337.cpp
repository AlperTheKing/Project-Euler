#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr int kDefaultN = 20'000'000;
constexpr u32 kMod = 100'000'000U;

struct Options {
    int n = kDefaultN;
    bool run_checks = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) return false;

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) return false;

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') return false;
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) return false;
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) return false;
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) return false;
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) return false;
    if (parsed > static_cast<u64>(std::numeric_limits<int>::max())) return false;
    value = static_cast<int>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checks") {
            options.run_checks = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        unsigned parsed_threads = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_threads)) {
            options.requested_threads = parsed_threads;
            continue;
        }

        int parsed_n = 0;
        if (parse_int_after_prefix(arg, "--n=", parsed_n)) {
            options.n = parsed_n;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n < 1) {
        std::cerr << "N must be positive.\n";
        return false;
    }
    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) return 1U;

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) threads = 1U;
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<int> primes;
    if (limit < 2) return primes;

    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(limit + 1), 0U);
    primes.reserve(static_cast<std::size_t>(limit / 10));

    for (int i = 2; i <= limit; ++i) {
        if (is_composite[static_cast<std::size_t>(i)] == 0U) {
            primes.push_back(i);
            if (static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(limit)) {
                for (u64 j = static_cast<u64>(i) * static_cast<u64>(i);
                     j <= static_cast<u64>(limit);
                     j += static_cast<u64>(i)) {
                    is_composite[static_cast<std::size_t>(j)] = 1U;
                }
            }
        }
    }

    return primes;
}

std::vector<u32> compute_totients_parallel(const int n, const unsigned thread_count) {
    std::vector<u32> phi(static_cast<std::size_t>(n + 1), 0U);
    if (n >= 1) phi[1] = 1U;

    const std::vector<int> primes = sieve_primes(n);
    if (n < 2) return phi;

    constexpr int kBlockSize = 1 << 20;
    const int block_count = (n + kBlockSize - 1) / kBlockSize;
    std::atomic<int> next_block{0};

    const auto process_block = [&](const int block_id) {
        const int left = block_id * kBlockSize + 1;
        const int right = std::min(n, (block_id + 1) * kBlockSize);

        for (int x = left; x <= right; ++x) {
            phi[static_cast<std::size_t>(x)] = static_cast<u32>(x);
        }

        for (int p : primes) {
            if (p > right) break;
            const int start = ((left + p - 1) / p) * p;
            for (int x = start; x <= right; x += p) {
                u32& cur = phi[static_cast<std::size_t>(x)];
                cur -= cur / static_cast<u32>(p);
            }
        }
    };

    if (thread_count <= 1U || block_count <= 1) {
        for (int block = 0; block < block_count; ++block) {
            process_block(block);
        }
        return phi;
    }

    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (unsigned t = 0; t < thread_count; ++t) {
        workers.emplace_back([&]() {
            while (true) {
                const int block_id = next_block.fetch_add(1, std::memory_order_relaxed);
                if (block_id >= block_count) break;
                process_block(block_id);
            }
        });
    }
    for (std::thread& th : workers) th.join();

    return phi;
}

class FenwickRangeAddPointQueryMod {
  public:
    explicit FenwickRangeAddPointQueryMod(const int n)
        : n_(n), bit_(static_cast<std::size_t>(n + 1), 0U) {}

    void range_add(int left, int right, const u32 delta) {
        if (left < 1) left = 1;
        if (right > n_) right = n_;
        if (left > right || delta == 0U) return;

        add_point(left, delta);
        if (right + 1 <= n_) {
            const u32 neg = (delta == 0U) ? 0U : (kMod - delta);
            add_point(right + 1, neg);
        }
    }

    u32 point_query(const int index) const {
        u32 out = 0U;
        int i = index;
        while (i > 0) {
            out += bit_[static_cast<std::size_t>(i)];
            if (out >= kMod) out -= kMod;
            i -= i & -i;
        }
        return out;
    }

  private:
    void add_point(int index, const u32 delta) {
        int i = index;
        while (i <= n_) {
            u32 next = bit_[static_cast<std::size_t>(i)] + delta;
            if (next >= kMod) next -= kMod;
            bit_[static_cast<std::size_t>(i)] = next;
            i += i & -i;
        }
    }

    int n_ = 0;
    std::vector<u32> bit_;
};

u32 count_sequences_mod(const int n, const std::vector<u32>& phi) {
    if (n < 6) return 0U;

    FenwickRangeAddPointQueryMod fenwick(n);
    u32 total = 1U;  // Sequence [6].
    fenwick.range_add(static_cast<int>(phi[6]) + 1, 5, 1U);

    for (int x = 7; x <= n; ++x) {
        const u32 dp_x = fenwick.point_query(static_cast<int>(phi[static_cast<std::size_t>(x)]));
        total += dp_x;
        if (total >= kMod) total -= kMod;
        fenwick.range_add(static_cast<int>(phi[static_cast<std::size_t>(x)]) + 1, x - 1, dp_x);
    }

    return total;
}

u64 count_sequences_exact_bruteforce(const int n, const std::vector<u32>& phi) {
    if (n < 6) return 0ULL;

    std::vector<u64> dp(static_cast<std::size_t>(n + 1), 0ULL);
    dp[6] = 1ULL;
    u64 total = 1ULL;

    for (int x = 7; x <= n; ++x) {
        const int phix = static_cast<int>(phi[static_cast<std::size_t>(x)]);
        u64 sum = 0ULL;
        for (int i = 6; i < x; ++i) {
            if (phi[static_cast<std::size_t>(i)] < static_cast<u32>(phix) && phix < i) {
                sum += dp[static_cast<std::size_t>(i)];
            }
        }
        dp[static_cast<std::size_t>(x)] = sum;
        total += sum;
    }

    return total;
}

bool validate_totients(const std::vector<u32>& phi) {
    if (phi.size() < 11U) return false;
    const std::vector<u32> expected = {0U, 1U, 1U, 2U, 2U, 4U, 2U, 6U, 4U, 6U, 4U};
    for (std::size_t i = 1; i <= 10U; ++i) {
        if (phi[i] != expected[i]) {
            std::cerr << "Totient validation failed at n=" << i << ": got " << phi[i]
                      << ", expected " << expected[i] << '\n';
            return false;
        }
    }
    return true;
}

bool run_validations(const bool allow_multithreading, const unsigned requested_threads) {
    {
        const std::vector<u32> phi_100 = compute_totients_parallel(100, 1U);
        if (!validate_totients(phi_100)) return false;

        const u64 s10 = count_sequences_exact_bruteforce(10, phi_100);
        if (s10 != 4ULL) {
            std::cerr << "Validation failed: S(10) = " << s10 << ", expected 4\n";
            return false;
        }

        const u64 s100 = count_sequences_exact_bruteforce(100, phi_100);
        if (s100 != 482'073'668ULL) {
            std::cerr << "Validation failed: S(100) = " << s100
                      << ", expected 482073668\n";
            return false;
        }

        const u32 s100_mod = count_sequences_mod(100, phi_100);
        if (s100_mod != static_cast<u32>(s100 % static_cast<u64>(kMod))) {
            std::cerr << "Validation failed: modular solver mismatch at N=100\n";
            return false;
        }
    }

    const unsigned threads_for_checks = choose_thread_count(
        allow_multithreading, requested_threads, std::max<std::size_t>(2ULL, 8ULL));

    {
        const std::vector<u32> phi_10k = compute_totients_parallel(10'000, threads_for_checks);
        const u32 s10k = count_sequences_mod(10'000, phi_10k);
        if (s10k != 73'808'307U) {
            std::cerr << "Validation failed: S(10000) mod 10^8 = " << s10k
                      << ", expected 73808307\n";
            return false;
        }
    }

    if (threads_for_checks > 1U) {
        const int test_n = 300'000;
        const std::vector<u32> phi_single = compute_totients_parallel(test_n, 1U);
        const std::vector<u32> phi_multi = compute_totients_parallel(test_n, threads_for_checks);
        if (phi_single != phi_multi) {
            std::cerr << "Validation failed: single-thread and multi-thread totients differ\n";
            return false;
        }

        const u32 s_single = count_sequences_mod(test_n, phi_single);
        const u32 s_multi = count_sequences_mod(test_n, phi_multi);
        if (s_single != s_multi) {
            std::cerr << "Validation failed: thread consistency mismatch in DP stage\n";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) return 1;

    if (options.run_checks && !run_validations(options.allow_multithreading, options.requested_threads)) {
        return 1;
    }

    const int n = options.n;
    const std::size_t totient_work_units =
        static_cast<std::size_t>((n + (1 << 20) - 1) / (1 << 20));
    const unsigned threads_for_totients = choose_thread_count(
        options.allow_multithreading, options.requested_threads, std::max<std::size_t>(1ULL, totient_work_units));

    const std::vector<u32> phi = compute_totients_parallel(n, threads_for_totients);
    const u32 answer = count_sequences_mod(n, phi);
    std::cout << answer << '\n';
    return 0;
}
