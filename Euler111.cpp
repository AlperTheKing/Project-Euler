#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int digits = 10;
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
        if (parse_int_after_prefix(arg, "--digits=", options.digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.digits >= 2;
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

bool is_prime(u64 n) {
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

    for (u64 a : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL}) {
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

void choose_positions(const int n,
                      const int k,
                      const int start,
                      std::vector<int>& current,
                      std::vector<std::vector<int>>& all) {
    if (static_cast<int>(current.size()) == k) {
        all.push_back(current);
        return;
    }

    for (int i = start; i <= n - (k - static_cast<int>(current.size())); ++i) {
        current.push_back(i);
        choose_positions(n, k, i + 1, current, all);
        current.pop_back();
    }
}

void assign_replacements(const int idx,
                         const std::vector<int>& positions,
                         const int repeated_digit,
                         std::vector<int>& digits,
                         std::vector<u64>& primes) {
    if (idx == static_cast<int>(positions.size())) {
        if (digits[0] == 0) {
            return;
        }
        const int last = digits.back();
        if ((last % 2) == 0 || last == 5) {
            return;
        }

        u64 value = 0;
        for (int d : digits) {
            value = value * 10ULL + static_cast<u64>(d);
        }
        if (is_prime(value)) {
            primes.push_back(value);
        }
        return;
    }

    for (int d = 0; d <= 9; ++d) {
        if (d == repeated_digit) {
            continue;
        }
        digits[static_cast<std::size_t>(positions[static_cast<std::size_t>(idx)])] = d;
        assign_replacements(idx + 1, positions, repeated_digit, digits, primes);
    }
}

u64 sum_for_digit(const int n, const int repeated_digit) {
    for (int replacements = 1; replacements <= n; ++replacements) {
        std::vector<std::vector<int>> position_sets;
        std::vector<int> current;
        choose_positions(n, replacements, 0, current, position_sets);

        std::vector<u64> primes;
        std::vector<int> digits(static_cast<std::size_t>(n), repeated_digit);

        for (const std::vector<int>& positions : position_sets) {
            assign_replacements(0, positions, repeated_digit, digits, primes);
            for (int pos : positions) {
                digits[static_cast<std::size_t>(pos)] = repeated_digit;
            }
        }

        if (!primes.empty()) {
            std::sort(primes.begin(), primes.end());
            primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
            u64 sum = 0;
            for (u64 p : primes) {
                sum += p;
            }
            return sum;
        }
    }

    return 0;
}

u64 solve(const int digits) {
    u64 total = 0;
    for (int d = 0; d <= 9; ++d) {
        total += sum_for_digit(digits, d);
    }
    return total;
}

bool run_checkpoints() {
    if (solve(4) != 273700ULL) {
        std::cerr << "Checkpoint failed for n=4" << '\n';
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

    std::cout << solve(options.digits) << '\n';
    return 0;
}
