#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <algorithm>
#include <functional>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int target = 4000000;
    bool run_checkpoints = true;
};

struct SearchState {
    u128 best_n = static_cast<u128>(~0ULL);
    int need_divisors = 0;
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
        if (parse_int_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.target >= 1;
}

void dfs_min_n(const std::array<int, 16>& primes,
               int idx,
               int max_exp,
               u128 current_n,
               int divisor_count,
               SearchState& state) {
    if (divisor_count >= state.need_divisors) {
        if (current_n < state.best_n) {
            state.best_n = current_n;
        }
        return;
    }

    if (idx >= static_cast<int>(primes.size())) {
        return;
    }

    u128 n = current_n;
    for (int exp = 1; exp <= max_exp; ++exp) {
        n *= static_cast<u128>(primes[static_cast<std::size_t>(idx)]);
        if (n >= state.best_n) {
            break;
        }

        const int new_divisor_count = divisor_count * (2 * exp + 1);
        dfs_min_n(primes, idx + 1, exp, n, new_divisor_count, state);
    }
}

u128 solve_u128(const int target) {
    static constexpr std::array<int, 16> primes = {
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37, 41, 43, 47, 53,
    };

    SearchState state;
    state.need_divisors = 2 * target + 1;
    dfs_min_n(primes, 0, 63, 1, 1, state);
    return state.best_n;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool run_checkpoints() {
    if (to_string_u128(solve_u128(4)) != "6") {
        std::cerr << "Checkpoint failed for target=4" << '\n';
        return false;
    }
    if (to_string_u128(solve_u128(1000)) != "180180") {
        std::cerr << "Checkpoint failed for target=1000" << '\n';
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

    std::cout << to_string_u128(solve_u128(options.target)) << '\n';
    return 0;
}
