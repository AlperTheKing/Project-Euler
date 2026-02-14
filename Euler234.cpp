#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
    u64 limit = 999966663333ULL;
    bool run_checkpoints = true;
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
        parsed = parsed * 10 + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 4;
}

std::vector<int> primes_up_to(const int n) {
    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;

    for (int i = 2; i <= n; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
        for (int p : primes) {
            const long long v = 1LL * i * p;
            if (v > n) {
                break;
            }
            is_composite[static_cast<std::size_t>(v)] = 1;
            if (i % p == 0) {
                break;
            }
        }
    }

    return primes;
}

i128 sum_multiples_in_range(const u64 d, const u64 left, const u64 right) {
    if (left > right) {
        return 0;
    }
    const u64 first = ((left + d - 1) / d) * d;
    const u64 last = (right / d) * d;
    if (first > last) {
        return 0;
    }

    const u64 count = (last - first) / d + 1;
    return (static_cast<i128>(first) + static_cast<i128>(last)) * static_cast<i128>(count) / 2;
}

std::string i128_to_string(i128 value) {
    if (value == 0) {
        return "0";
    }

    bool neg = value < 0;
    if (neg) {
        value = -value;
    }

    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (neg) {
        out.push_back('-');
    }
    std::reverse(out.begin(), out.end());
    return out;
}

i128 solve_i128(const u64 limit) {
    const int sqrt_limit = static_cast<int>(std::sqrt(static_cast<long double>(limit))) + 10;
    std::vector<int> primes = primes_up_to(sqrt_limit + 100);

    // Ensure one prime above sqrt(limit) exists.
    int candidate = sqrt_limit + 101;
    while (true) {
        bool is_prime = true;
        for (int p : primes) {
            if (1LL * p * p > candidate) {
                break;
            }
            if (candidate % p == 0) {
                is_prime = false;
                break;
            }
        }
        if (is_prime) {
            primes.push_back(candidate);
            break;
        }
        ++candidate;
    }

    i128 total = 0;
    for (std::size_t i = 0; i + 1 < primes.size(); ++i) {
        const u64 p = static_cast<u64>(primes[i]);
        const u64 q = static_cast<u64>(primes[i + 1]);

        const u64 p2 = p * p;
        if (p2 > limit) {
            break;
        }

        const u64 q2 = q * q;
        const u64 left = p2 + 1;
        const u64 right = std::min(limit, q2 - 1);
        if (left > right) {
            continue;
        }

        const i128 sum_p = sum_multiples_in_range(p, left, right);
        const i128 sum_q = sum_multiples_in_range(q, left, right);
        const i128 sum_pq = sum_multiples_in_range(p * q, left, right);
        total += sum_p + sum_q - 2 * sum_pq;
    }

    return total;
}

u64 solve_u64(const u64 limit) {
    return static_cast<u64>(solve_i128(limit));
}

bool run_checkpoints() {
    if (solve_u64(15ULL) != 30ULL) {
        std::cerr << "Checkpoint failed for limit 15" << '\n';
        return false;
    }
    if (solve_u64(1000ULL) != 34825ULL) {
        std::cerr << "Checkpoint failed for limit 1000" << '\n';
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

    std::cout << i128_to_string(solve_i128(options.limit)) << '\n';
    return 0;
}
