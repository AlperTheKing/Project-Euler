#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

struct Options {
    int family_size = 8;
    int sieve_limit = 2000000;
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
        if (parse_int_after_prefix(arg, "--family-size=", options.family_size)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--sieve-limit=", options.sieve_limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.family_size >= 2 && options.sieve_limit >= 100;
}

std::vector<bool> sieve(const int limit) {
    std::vector<bool> prime(static_cast<std::size_t>(limit + 1), true);
    prime[0] = false;
    prime[1] = false;

    for (int p = 2; p <= limit / p; ++p) {
        if (!prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            prime[static_cast<std::size_t>(q)] = false;
        }
    }

    return prime;
}

bool is_prime_trial(std::int64_t n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (std::int64_t p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

bool is_prime_fast(const int n, const std::vector<bool>& sieve_table) {
    if (n < static_cast<int>(sieve_table.size())) {
        return sieve_table[static_cast<std::size_t>(n)];
    }
    return is_prime_trial(n);
}

int solve(const int family_size, const int sieve_limit) {
    const std::vector<bool> sieve_table = sieve(sieve_limit);

    for (int p = 11; p <= sieve_limit; ++p) {
        if (!sieve_table[static_cast<std::size_t>(p)]) {
            continue;
        }

        const std::string base = std::to_string(p);

        for (char target_digit = '0'; target_digit <= '9'; ++target_digit) {
            std::vector<int> positions;
            for (int i = 0; i < static_cast<int>(base.size()); ++i) {
                if (base[static_cast<std::size_t>(i)] == target_digit) {
                    positions.push_back(i);
                }
            }
            if (positions.empty()) {
                continue;
            }

            const int subsets = 1 << static_cast<int>(positions.size());
            for (int mask = 1; mask < subsets; ++mask) {
                int family_count = 0;
                int smallest_family_prime = std::numeric_limits<int>::max();

                for (char replacement = '0'; replacement <= '9'; ++replacement) {
                    std::string candidate = base;
                    for (int bit = 0; bit < static_cast<int>(positions.size()); ++bit) {
                        if ((mask >> bit) & 1) {
                            candidate[static_cast<std::size_t>(positions[static_cast<std::size_t>(bit)])] =
                                replacement;
                        }
                    }

                    if (candidate[0] == '0') {
                        continue;
                    }

                    const int value = std::stoi(candidate);
                    if (is_prime_fast(value, sieve_table)) {
                        ++family_count;
                        if (value < smallest_family_prime) {
                            smallest_family_prime = value;
                        }
                    }
                }

                if (family_count >= family_size && smallest_family_prime == p) {
                    return p;
                }
            }
        }
    }

    return 0;
}

bool run_checkpoints() {
    if (solve(6, 100000) != 13) {
        std::cerr << "Checkpoint failed for family size 6" << '\n';
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

    std::cout << solve(options.family_size, options.sieve_limit) << '\n';
    return 0;
}
