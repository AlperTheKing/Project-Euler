#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int from = 10000000;
    int to = 20000000;
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
        if (parse_int_after_prefix(arg, "--from=", options.from) ||
            parse_int_after_prefix(arg, "--to=", options.to)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.from >= 2 && options.to > options.from;
}

std::array<int, 10> digit_masks() {
    // Segment order: top, middle, bottom, top-left, top-right, bottom-left, bottom-right.
    constexpr int T = 1 << 0;
    constexpr int M = 1 << 1;
    constexpr int B = 1 << 2;
    constexpr int TL = 1 << 3;
    constexpr int TR = 1 << 4;
    constexpr int BL = 1 << 5;
    constexpr int BR = 1 << 6;

    return {
        T | B | TL | TR | BL | BR,          // 0
        TR | BR,                            // 1
        T | M | B | TR | BL,                // 2
        T | M | B | TR | BR,                // 3
        M | TL | TR | BR,                   // 4
        T | M | B | TL | BR,                // 5
        T | M | B | TL | BL | BR,           // 6
        T | TL | TR | BR,                   // 7 (4 segments in this display style)
        T | M | B | TL | TR | BL | BR,      // 8
        T | M | B | TL | TR | BR            // 9
    };
}

int digit_sum(int n) {
    int s = 0;
    while (n > 0) {
        s += n % 10;
        n /= 10;
    }
    return s;
}

std::vector<int> digital_root_chain(int n) {
    std::vector<int> chain;
    while (true) {
        chain.push_back(n);
        if (n < 10) {
            break;
        }
        n = digit_sum(n);
    }
    return chain;
}

int segment_count_number(int n, const std::array<int, 10>& masks) {
    int total = 0;
    while (n > 0) {
        total += __builtin_popcount(static_cast<unsigned>(masks[static_cast<std::size_t>(n % 10)]));
        n /= 10;
    }
    return total;
}

int transition_cost(int from, int to, const std::array<int, 10>& masks) {
    int cost = 0;
    while (from > 0 || to > 0) {
        const int a = (from > 0) ? masks[static_cast<std::size_t>(from % 10)] : 0;
        const int b = (to > 0) ? masks[static_cast<std::size_t>(to % 10)] : 0;
        cost += __builtin_popcount(static_cast<unsigned>(a ^ b));
        from /= 10;
        to /= 10;
    }
    return cost;
}

int sam_cost(const std::vector<int>& chain, const std::array<int, 10>& masks) {
    int cost = 0;
    for (int value : chain) {
        cost += 2 * segment_count_number(value, masks);
    }
    return cost;
}

int max_cost(const std::vector<int>& chain, const std::array<int, 10>& masks) {
    int cost = 0;
    cost += transition_cost(0, chain.front(), masks);
    for (std::size_t i = 1; i < chain.size(); ++i) {
        cost += transition_cost(chain[i - 1], chain[i], masks);
    }
    cost += transition_cost(chain.back(), 0, masks);
    return cost;
}

i64 solve(const int from, const int to) {
    const auto masks = digit_masks();

    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(to), 1U);
    is_prime[0] = 0U;
    if (to > 1) {
        is_prime[1] = 0U;
    }
    for (int p = 2; static_cast<long long>(p) * p < to; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (int q = p * p; q < to; q += p) {
            is_prime[static_cast<std::size_t>(q)] = 0U;
        }
    }

    i64 total_diff = 0;
    for (int p = from; p < to; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        const auto chain = digital_root_chain(p);
        total_diff += static_cast<i64>(sam_cost(chain, masks) - max_cost(chain, masks));
    }
    return total_diff;
}

bool run_checkpoints() {
    const auto masks = digit_masks();
    {
        const auto chain = digital_root_chain(137);
        const int sam = sam_cost(chain, masks);
        const int mx = max_cost(chain, masks);
        if (sam != 40 || mx != 30) {
            std::cerr << "Checkpoint failed for 137 sample transitions" << '\n';
            return false;
        }
    }
    {
        const auto chain = digital_root_chain(7);
        if (sam_cost(chain, masks) != max_cost(chain, masks)) {
            std::cerr << "Checkpoint failed for single-digit transition equivalence" << '\n';
            return false;
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
    std::cout << solve(options.from, options.to) << '\n';
    return 0;
}
