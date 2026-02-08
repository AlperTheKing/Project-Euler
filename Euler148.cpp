#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 rows = 1000000000ULL;
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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--rows=", options.rows)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

u128 count_not_divisible_by_seven(const u64 rows) {
    if (rows == 0) {
        return 0;
    }

    std::vector<int> digits;
    u64 x = rows;
    while (x > 0) {
        digits.push_back(static_cast<int>(x % 7ULL));
        x /= 7ULL;
    }
    std::reverse(digits.begin(), digits.end());

    std::vector<u128> pow28(static_cast<std::size_t>(digits.size() + 1), 1);
    for (std::size_t i = 1; i < pow28.size(); ++i) {
        pow28[i] = pow28[i - 1] * 28U;
    }

    u128 answer = 0;
    u128 prefix_product = 1;

    for (std::size_t i = 0; i < digits.size(); ++i) {
        const int d = digits[i];
        const std::size_t rem = digits.size() - i - 1;

        for (int xdigit = 0; xdigit < d; ++xdigit) {
            answer += prefix_product * static_cast<u128>(xdigit + 1) * pow28[rem];
        }
        prefix_product *= static_cast<u128>(d + 1);
    }

    return answer;
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
    if (to_string_u128(count_not_divisible_by_seven(100ULL)) != "2361") {
        std::cerr << "Checkpoint failed for first 100 rows" << '\n';
        return false;
    }
    if (to_string_u128(count_not_divisible_by_seven(7ULL)) != "28") {
        std::cerr << "Checkpoint failed for first 7 rows" << '\n';
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

    std::cout << to_string_u128(count_not_divisible_by_seven(options.rows)) << '\n';
    return 0;
}
