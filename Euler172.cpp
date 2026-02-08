#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u128 = unsigned __int128;

struct Options {
    int length = 18;
    int max_repeats = 3;
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
        if (parse_int_after_prefix(arg, "--length=", options.length) ||
            parse_int_after_prefix(arg, "--max-repeats=", options.max_repeats)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.length >= 1 && options.max_repeats >= 1;
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

u128 solve(const int length, const int max_repeats) {
    std::vector<u128> fact(static_cast<std::size_t>(length + 1), 1);
    for (int i = 1; i <= length; ++i) {
        fact[static_cast<std::size_t>(i)] = fact[static_cast<std::size_t>(i - 1)] * static_cast<u128>(i);
    }

    std::array<int, 10> counts{};
    u128 total = 0;

    const auto dfs = [&](auto&& self, int digit, int remaining, u128 denom) -> void {
        if (digit == 10) {
            if (remaining != 0) {
                return;
            }

            const u128 all_permutations = fact[static_cast<std::size_t>(length)] / denom;
            u128 leading_zero = 0;
            if (counts[0] > 0) {
                const u128 denom2 = (denom / fact[static_cast<std::size_t>(counts[0])]) *
                                    fact[static_cast<std::size_t>(counts[0] - 1)];
                leading_zero = fact[static_cast<std::size_t>(length - 1)] / denom2;
            }

            total += all_permutations - leading_zero;
            return;
        }

        const int max_take = std::min(max_repeats, remaining);
        for (int c = 0; c <= max_take; ++c) {
            counts[static_cast<std::size_t>(digit)] = c;
            self(self, digit + 1, remaining - c, denom * fact[static_cast<std::size_t>(c)]);
        }
    };

    dfs(dfs, 0, length, 1);
    return total;
}

bool run_checkpoints() {
    if (to_string_u128(solve(4, 3)) != "8991") {
        std::cerr << "Checkpoint failed for length=4" << '\n';
        return false;
    }
    if (to_string_u128(solve(3, 3)) != "900") {
        std::cerr << "Checkpoint failed for length=3" << '\n';
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

    std::cout << to_string_u128(solve(options.length, options.max_repeats)) << '\n';
    return 0;
}
