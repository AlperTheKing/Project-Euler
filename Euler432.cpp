#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kN = 510'510ULL;
constexpr u64 kCheckpointM = 1'000'000ULL;
constexpr u64 kTargetM = 100'000'000'000ULL;
constexpr u64 kExpectedCheckpoint = 45'480'596'821'125'120ULL;
constexpr u64 kLastDigitsMod = 1'000'000'000ULL;
constexpr u64 kTotientSieveLimit = 12'000'000ULL;

struct U64Hash {
    std::size_t operator()(u64 value) const noexcept {
        value += 0x9e3779b97f4a7c15ULL;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        value ^= (value >> 31);
        return static_cast<std::size_t>(value);
    }
};

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(digits.begin(), digits.end());
    return digits;
}

u128 triangular_u128(u64 value) {
    return static_cast<u128>(value) * (value + 1) / 2;
}

u64 euler_phi_u64(u64 value) {
    if (value == 0) {
        return 0;
    }

    u64 result = value;
    u64 n = value;

    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        while (n % p == 0) {
            n /= p;
        }
        result -= result / p;
    }
    if (n > 1) {
        result -= result / n;
    }
    return result;
}

std::vector<u64> distinct_prime_factors(u64 n) {
    std::vector<u64> factors;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        factors.push_back(p);
        while (n % p == 0) {
            n /= p;
        }
    }
    if (n > 1) {
        factors.push_back(n);
    }
    return factors;
}

void generate_smooth_numbers_dfs(const std::vector<u64>& primes,
                                 std::size_t index,
                                 u64 current,
                                 u64 limit,
                                 std::vector<u64>& out) {
    if (index == primes.size()) {
        out.push_back(current);
        return;
    }

    u64 value = current;
    const u64 prime = primes[index];
    while (value <= limit) {
        generate_smooth_numbers_dfs(primes, index + 1, value, limit, out);
        if (value > limit / prime) {
            break;
        }
        value *= prime;
    }
}

std::vector<u64> generate_smooth_numbers(const std::vector<u64>& primes, u64 limit) {
    std::vector<u64> smooth;
    smooth.reserve(200'000);
    generate_smooth_numbers_dfs(primes, 0, 1, limit, smooth);
    return smooth;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 100'000) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    threads = std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
    return threads;
}

class TotientSummatory {
public:
    explicit TotientSummatory(u64 sieve_limit)
        : sieve_limit_(sieve_limit),
          phi_(static_cast<std::size_t>(sieve_limit_) + 1, 0),
          prefix_(static_cast<std::size_t>(sieve_limit_) + 1, 0) {
        build_totients();
        large_cache_.reserve(1 << 20);
    }

    u128 sum_phi(u64 n) {
        if (n <= sieve_limit_) {
            return prefix_[static_cast<std::size_t>(n)];
        }

        const auto it = large_cache_.find(n);
        if (it != large_cache_.end()) {
            return it->second;
        }

        u128 result = triangular_u128(n);
        for (u64 left = 2, right = 0; left <= n; left = right + 1) {
            const u64 quotient = n / left;
            right = n / quotient;
            result -= static_cast<u128>(right - left + 1) * sum_phi(quotient);
        }

        large_cache_.emplace(n, result);
        return result;
    }

private:
    u64 sieve_limit_ = 0;
    std::vector<u32> phi_;
    std::vector<u64> prefix_;
    std::unordered_map<u64, u128, U64Hash> large_cache_;

    void build_totients() {
        std::vector<int> primes;
        primes.reserve(static_cast<std::size_t>(sieve_limit_ / 10));
        std::vector<char> is_composite(static_cast<std::size_t>(sieve_limit_) + 1, 0);

        if (sieve_limit_ >= 1) {
            phi_[1] = 1;
        }

        for (u64 i = 2; i <= sieve_limit_; ++i) {
            if (!is_composite[static_cast<std::size_t>(i)]) {
                primes.push_back(static_cast<int>(i));
                phi_[static_cast<std::size_t>(i)] = static_cast<u32>(i - 1);
            }

            for (const int prime_int : primes) {
                const u64 prime = static_cast<u64>(prime_int);
                const u64 composite = i * prime;
                if (composite > sieve_limit_) {
                    break;
                }

                is_composite[static_cast<std::size_t>(composite)] = 1;
                if (i % prime == 0) {
                    phi_[static_cast<std::size_t>(composite)] =
                        static_cast<u32>(phi_[static_cast<std::size_t>(i)] * prime);
                    break;
                }

                phi_[static_cast<std::size_t>(composite)] = static_cast<u32>(
                    static_cast<u64>(phi_[static_cast<std::size_t>(i)]) * (prime - 1));
            }
        }

        for (u64 i = 1; i <= sieve_limit_; ++i) {
            prefix_[static_cast<std::size_t>(i)] =
                prefix_[static_cast<std::size_t>(i - 1)] +
                static_cast<u64>(phi_[static_cast<std::size_t>(i)]);
        }
    }
};

class Euler432Solver {
public:
    explicit Euler432Solver(u64 n)
        : n_(n),
          prime_factors_(distinct_prime_factors(n)),
          phi_n_(euler_phi_u64(n)),
          totient_summatory_(kTotientSieveLimit) {}

    u128 solve(u64 m, bool allow_multithreading, unsigned requested_threads = 0) {
        std::vector<u64> smooth = generate_smooth_numbers(prime_factors_, m);

        std::vector<u64> quotients;
        quotients.reserve(smooth.size());
        for (const u64 value : smooth) {
            quotients.push_back(m / value);
        }

        std::vector<u64> unique_quotients = quotients;
        std::sort(unique_quotients.begin(), unique_quotients.end());
        unique_quotients.erase(std::unique(unique_quotients.begin(), unique_quotients.end()),
                               unique_quotients.end());

        std::unordered_map<u64, u128, U64Hash> sum_phi_by_quotient;
        sum_phi_by_quotient.reserve(unique_quotients.size() * 2 + 16);
        for (const u64 q : unique_quotients) {
            sum_phi_by_quotient.emplace(q, totient_summatory_.sum_phi(q));
        }

        const u128 inner = sum_over_quotients(quotients,
                                              sum_phi_by_quotient,
                                              allow_multithreading,
                                              requested_threads);
        return static_cast<u128>(phi_n_) * inner;
    }

private:
    u64 n_ = 0;
    std::vector<u64> prime_factors_;
    u64 phi_n_ = 0;
    TotientSummatory totient_summatory_;

    static u128 sum_over_quotients(const std::vector<u64>& quotients,
                                   const std::unordered_map<u64, u128, U64Hash>& sum_phi_by_quotient,
                                   bool allow_multithreading,
                                   unsigned requested_threads) {
        const unsigned threads =
            choose_thread_count(allow_multithreading, requested_threads, quotients.size());

        if (threads == 1) {
            u128 total = 0;
            for (const u64 q : quotients) {
                total += sum_phi_by_quotient.at(q);
            }
            return total;
        }

        std::vector<u128> partial(threads, 0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                const std::size_t left = quotients.size() * tid / threads;
                const std::size_t right = quotients.size() * (tid + 1) / threads;

                u128 local = 0;
                for (std::size_t i = left; i < right; ++i) {
                    local += sum_phi_by_quotient.at(quotients[i]);
                }
                partial[tid] = local;
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        u128 total = 0;
        for (const u128 value : partial) {
            total += value;
        }
        return total;
    }
};

u128 brute_force_s(u64 n, u64 m) {
    u128 total = 0;
    for (u64 i = 1; i <= m; ++i) {
        total += euler_phi_u64(n * i);
    }
    return total;
}

bool parse_threads_argument(const std::string& arg, unsigned& threads) {
    constexpr const char* prefix = "--threads=";
    if (arg.rfind(prefix, 0) != 0) {
        return false;
    }

    const std::string value = arg.substr(std::char_traits<char>::length(prefix));
    if (value.empty()) {
        return false;
    }

    unsigned parsed = 0;
    for (const char c : value) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<unsigned>(c - '0');
    }

    threads = parsed;
    return true;
}

bool run_validation_checkpoints(Euler432Solver& solver,
                                bool allow_multithreading,
                                unsigned requested_threads) {
    const std::vector<u64> brute_small = {1, 2, 5, 10, 100, 1'000, 10'000};

    for (const u64 m : brute_small) {
        const u128 fast = solver.solve(m, false, 1);
        const u128 brute = brute_force_s(kN, m);
        if (fast != brute) {
            std::cerr << "Brute-force checkpoint failed for m=" << m
                      << ": fast=" << to_string_u128(fast)
                      << ", brute=" << to_string_u128(brute) << '\n';
            return false;
        }
        std::cout << "Brute checkpoint S(510510, " << m << ") = "
                  << to_string_u128(fast) << " (ok)\n";
    }

    const u128 given = solver.solve(kCheckpointM, allow_multithreading, requested_threads);
    if (given != static_cast<u128>(kExpectedCheckpoint)) {
        std::cerr << "Given checkpoint failed for S(510510, 10^6): expected "
                  << kExpectedCheckpoint << ", got " << to_string_u128(given) << '\n';
        return false;
    }

    std::cout << "Given checkpoint S(510510, 10^6) = " << to_string_u128(given)
              << " (ok)\n";
    return true;
}

} // namespace

int main(int argc, char** argv) {
    bool allow_multithreading = true;
    unsigned requested_threads = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--single-thread") {
            allow_multithreading = false;
            continue;
        }
        if (!parse_threads_argument(arg, requested_threads)) {
            std::cerr << "Unknown argument: " << arg << '\n';
            std::cerr << "Usage: ./Euler432 [--single-thread] [--threads=N]\n";
            return 1;
        }
    }

    Euler432Solver solver(kN);
    if (!run_validation_checkpoints(solver, allow_multithreading, requested_threads)) {
        return 1;
    }

    const u128 answer = solver.solve(kTargetM, allow_multithreading, requested_threads);
    const u64 last_nine_digits = static_cast<u64>(answer % kLastDigitsMod);

    std::cout << "S(510510, 10^11) = " << to_string_u128(answer) << '\n';
    std::cout << "Last 9 digits: " << std::setw(9) << std::setfill('0') << last_nine_digits
              << '\n';
    return 0;
}
