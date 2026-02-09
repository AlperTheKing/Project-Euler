#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kMod = 1'000'000'007U;

struct Options {
    int n = 10'000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n <= 0) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

int estimate_limit_for_prime_count(const int n) {
    if (n < 6) {
        return 32;
    }
    const double dn = static_cast<double>(n);
    const double estimate = dn * (std::log(dn) + std::log(std::log(dn)));
    return static_cast<int>(estimate) + 64;
}

std::vector<bool> sieve(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    is_prime[0] = false;
    is_prime[1] = false;
    for (int p = 2; static_cast<int64_t>(p) * p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int x = p * p; x <= limit; x += p) {
            is_prime[static_cast<std::size_t>(x)] = false;
        }
    }
    return is_prime;
}

std::uint8_t digital_root(const int x) {
    return static_cast<std::uint8_t>(1 + (x - 1) % 9);
}

void build_sequences(const int n, std::vector<std::uint8_t>& prime_digits,
                     std::vector<std::uint8_t>& composite_digits) {
    int limit = std::max(64, estimate_limit_for_prime_count(n));
    while (true) {
        const std::vector<bool> is_prime = sieve(limit);
        prime_digits.clear();
        composite_digits.clear();
        prime_digits.reserve(static_cast<std::size_t>(n));
        composite_digits.reserve(static_cast<std::size_t>(n));

        for (int x = 2; x <= limit; ++x) {
            if (is_prime[static_cast<std::size_t>(x)]) {
                if (static_cast<int>(prime_digits.size()) < n) {
                    prime_digits.push_back(digital_root(x));
                }
            } else {
                if (static_cast<int>(composite_digits.size()) < n) {
                    composite_digits.push_back(digital_root(x));
                }
            }
            if (static_cast<int>(prime_digits.size()) == n &&
                static_cast<int>(composite_digits.size()) == n) {
                return;
            }
        }
        limit *= 2;
    }
}

inline std::size_t idx(const int i, const int j, const int width) {
    return static_cast<std::size_t>(i) * static_cast<std::size_t>(width) +
           static_cast<std::size_t>(j);
}

u32 compute_scs_mod(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b,
                    std::string* out_string = nullptr) {
    const int n = static_cast<int>(a.size());
    const int width = n + 1;
    std::vector<std::uint16_t> dp(static_cast<std::size_t>(width) * static_cast<std::size_t>(width),
                                  0U);

    for (int i = n; i >= 0; --i) {
        for (int j = n; j >= 0; --j) {
            const std::size_t at = idx(i, j, width);
            if (i == n) {
                dp[at] = static_cast<std::uint16_t>(n - j);
                continue;
            }
            if (j == n) {
                dp[at] = static_cast<std::uint16_t>(n - i);
                continue;
            }
            if (a[static_cast<std::size_t>(i)] == b[static_cast<std::size_t>(j)]) {
                dp[at] = static_cast<std::uint16_t>(1 + dp[idx(i + 1, j + 1, width)]);
            } else {
                const std::uint16_t da = dp[idx(i + 1, j, width)];
                const std::uint16_t db = dp[idx(i, j + 1, width)];
                dp[at] = static_cast<std::uint16_t>(1 + (da < db ? da : db));
            }
        }
    }

    int i = 0;
    int j = 0;
    u32 result_mod = 0U;
    if (out_string != nullptr) {
        out_string->clear();
        out_string->reserve(static_cast<std::size_t>(dp[0]));
    }

    while (i < n || j < n) {
        std::uint8_t digit = 0U;
        if (i == n) {
            digit = b[static_cast<std::size_t>(j++)];
        } else if (j == n) {
            digit = a[static_cast<std::size_t>(i++)];
        } else if (a[static_cast<std::size_t>(i)] == b[static_cast<std::size_t>(j)]) {
            digit = a[static_cast<std::size_t>(i)];
            ++i;
            ++j;
        } else {
            const std::uint16_t da = dp[idx(i + 1, j, width)];
            const std::uint16_t db = dp[idx(i, j + 1, width)];
            if (da < db || (da == db && a[static_cast<std::size_t>(i)] < b[static_cast<std::size_t>(j)])) {
                digit = a[static_cast<std::size_t>(i++)];
            } else {
                digit = b[static_cast<std::size_t>(j++)];
            }
        }

        result_mod = static_cast<u32>((static_cast<u64>(result_mod) * 10ULL + digit) % kMod);
        if (out_string != nullptr) {
            out_string->push_back(static_cast<char>('0' + digit));
        }
    }

    return result_mod;
}

u32 solve_mod(const int n, std::string* out_string = nullptr) {
    std::vector<std::uint8_t> prime_digits;
    std::vector<std::uint8_t> composite_digits;
    build_sequences(n, prime_digits, composite_digits);
    return compute_scs_mod(prime_digits, composite_digits, out_string);
}

bool run_checkpoints() {
    std::string f10;
    solve_mod(10, &f10);
    if (f10 != "2357246891352679") {
        std::cerr << "Checkpoint failed: f(10)\n";
        return false;
    }
    if (solve_mod(100, nullptr) != 771'661'825U) {
        std::cerr << "Checkpoint failed: f(100) mod 1e9+7\n";
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
        return 1;
    }

    std::cout << solve_mod(options.n, nullptr) << '\n';
    return 0;
}
