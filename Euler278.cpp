#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int prime_limit = 5000;
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
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_int_after_prefix(arg, "--prime-limit=", options.prime_limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.prime_limit >= 3;
}

std::vector<int> primes_below(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit), true);
    if (limit > 0) {
        is_prime[0] = false;
    }
    if (limit > 1) {
        is_prime[1] = false;
    }

    for (int p = 2; p * p < limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int x = p * p; x < limit; x += p) {
            is_prime[static_cast<std::size_t>(x)] = false;
        }
    }

    std::vector<int> primes;
    for (int p = 2; p < limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            primes.push_back(p);
        }
    }
    return primes;
}

u64 frobenius_semiprime_triple(const u64 p, const u64 q, const u64 r) {
    return 2ULL * p * q * r - p * q - p * r - q * r;
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 solve_sum(const int prime_limit) {
    const std::vector<int> primes = primes_below(prime_limit);

    u128 answer = 0;
    const int n = static_cast<int>(primes.size());
    for (int i = 0; i < n; ++i) {
        const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(i)]);
        for (int j = i + 1; j < n; ++j) {
            const u64 q = static_cast<u64>(primes[static_cast<std::size_t>(j)]);
            for (int k = j + 1; k < n; ++k) {
                const u64 r = static_cast<u64>(primes[static_cast<std::size_t>(k)]);
                answer += static_cast<u128>(frobenius_semiprime_triple(p, q, r));
            }
        }
    }
    return answer;
}

u64 brute_frobenius_three(const u64 a, const u64 b, const u64 c) {
    const u64 cap = 4 * a * b * c;
    std::vector<bool> can(static_cast<std::size_t>(cap + 1), false);
    can[0] = true;

    for (u64 x = 0; x <= cap; ++x) {
        if (!can[static_cast<std::size_t>(x)]) {
            continue;
        }
        if (x + a <= cap) {
            can[static_cast<std::size_t>(x + a)] = true;
        }
        if (x + b <= cap) {
            can[static_cast<std::size_t>(x + b)] = true;
        }
        if (x + c <= cap) {
            can[static_cast<std::size_t>(x + c)] = true;
        }
    }

    u64 best = 0;
    for (u64 x = 1; x <= cap; ++x) {
        if (!can[static_cast<std::size_t>(x)]) {
            best = x;
        }
    }
    return best;
}

bool run_checkpoints() {
    if (frobenius_semiprime_triple(2, 3, 5) != 29ULL) {
        std::cerr << "Checkpoint failed for (2,3,5)" << '\n';
        return false;
    }
    if (frobenius_semiprime_triple(2, 7, 11) != 195ULL) {
        std::cerr << "Checkpoint failed for (2,7,11)" << '\n';
        return false;
    }

    // Verify formula against brute force on a small prime pool.
    const std::vector<int> small_primes = primes_below(20);
    for (std::size_t i = 0; i < small_primes.size(); ++i) {
        for (std::size_t j = i + 1; j < small_primes.size(); ++j) {
            for (std::size_t k = j + 1; k < small_primes.size(); ++k) {
                const u64 p = static_cast<u64>(small_primes[i]);
                const u64 q = static_cast<u64>(small_primes[j]);
                const u64 r = static_cast<u64>(small_primes[k]);
                const u64 a = p * q;
                const u64 b = p * r;
                const u64 c = q * r;
                const u64 brute = brute_frobenius_three(a, b, c);
                const u64 formula = frobenius_semiprime_triple(p, q, r);
                if (brute != formula) {
                    std::cerr << "Formula mismatch for (" << p << ',' << q << ',' << r
                              << "): brute=" << brute << ", formula=" << formula << '\n';
                    return false;
                }
            }
        }
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

    std::cout << to_string_u128(solve_sum(options.prime_limit)) << '\n';
    return 0;
}
