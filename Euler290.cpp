#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;

struct Options {
    int digits = 18;
    int multiplier = 137;
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
        if (parse_int_after_prefix(arg, "--digits=", options.digits) ||
            parse_int_after_prefix(arg, "--multiplier=", options.multiplier)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.digits >= 1 && options.multiplier >= 1;
}

int digit_sum(u64 x) {
    int s = 0;
    while (x > 0ULL) {
        s += static_cast<int>(x % 10ULL);
        x /= 10ULL;
    }
    return s;
}

u64 solve(const int digits, const int multiplier) {
    // key = carry * 1000 + (diff + offset)
    // diff = sum(processed digits of n) - sum(processed low digits of multiplier*n)
    static constexpr int kOffset = 500;

    std::unordered_map<int, u64> dp;
    dp.reserve(100000);
    dp[0 * 1000 + kOffset] = 1ULL;

    for (int pos = 0; pos < digits; ++pos) {
        std::unordered_map<int, u64> next;
        next.reserve(dp.size() * 8U);
        for (const auto& [key, ways] : dp) {
            const int carry = key / 1000;
            const int diff = (key % 1000) - kOffset;
            for (int d = 0; d <= 9; ++d) {
                const int t = multiplier * d + carry;
                const int out_digit = t % 10;
                const int next_carry = t / 10;
                const int next_diff = diff + d - out_digit;
                const int next_key = next_carry * 1000 + (next_diff + kOffset);
                next[next_key] += ways;
            }
        }
        dp.swap(next);
    }

    u64 answer = 0ULL;
    for (const auto& [key, ways] : dp) {
        const int carry = key / 1000;
        const int diff = (key % 1000) - kOffset;
        if (diff == digit_sum(static_cast<u64>(carry))) {
            answer += ways;
        }
    }
    return answer;
}

u64 brute_small(const int digits, const int multiplier) {
    u64 upper = 1ULL;
    for (int i = 0; i < digits; ++i) {
        upper *= 10ULL;
    }
    u64 count = 0ULL;
    for (u64 n = 0; n < upper; ++n) {
        if (digit_sum(n) == digit_sum(static_cast<u64>(multiplier) * n)) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(4, 137) != brute_small(4, 137)) {
        std::cerr << "Checkpoint failed for digits=4 multiplier=137" << '\n';
        return false;
    }
    if (solve(5, 137) != brute_small(5, 137)) {
        std::cerr << "Checkpoint failed for digits=5 multiplier=137" << '\n';
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
    std::cout << solve(options.digits, options.multiplier) << '\n';
    return 0;
}
