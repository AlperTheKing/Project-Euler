#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int n = 100;
    int choose = 50;
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
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--choose=", options.choose)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.n >= 1 && options.choose >= 0 && options.choose <= options.n;
}

i64 solve(const int n, const int choose) {
    std::vector<int> values;
    values.reserve(static_cast<std::size_t>(n));
    for (int i = 1; i <= n; ++i) {
        values.push_back(i * i);
    }

    int max_sum = 0;
    for (int i = n - choose; i < n; ++i) {
        if (i >= 0) {
            max_sum += values[static_cast<std::size_t>(i)];
        }
    }

    const int width = max_sum + 1;
    std::vector<std::uint8_t> dp(static_cast<std::size_t>((choose + 1) * width), 0);

    auto idx = [&](const int k, const int sum) -> std::size_t {
        return static_cast<std::size_t>(k) * static_cast<std::size_t>(width) +
               static_cast<std::size_t>(sum);
    };

    dp[idx(0, 0)] = 1;

    std::vector<int> low(choose + 1, 0);
    std::vector<int> high(choose + 1, -1);
    high[0] = 0;

    for (int i = 0; i < n; ++i) {
        const int x = values[static_cast<std::size_t>(i)];
        const int up = std::min(choose, i + 1);

        for (int k = up; k >= 1; --k) {
            if (high[k - 1] < 0) {
                continue;
            }

            const int start = low[k - 1] + x;
            const int end = high[k - 1] + x;

            for (int s = end; s >= start; --s) {
                const std::uint8_t add = dp[idx(k - 1, s - x)];
                if (add == 0) {
                    continue;
                }

                const std::size_t p = idx(k, s);
                const int sum = static_cast<int>(dp[p]) + static_cast<int>(add);
                dp[p] = static_cast<std::uint8_t>(sum >= 2 ? 2 : sum);
            }

            if (high[k] < 0) {
                low[k] = start;
                high[k] = end;
            } else {
                low[k] = std::min(low[k], start);
                high[k] = std::max(high[k], end);
            }
        }
    }

    i64 answer = 0;
    if (choose == 0) {
        return 0;
    }
    for (int s = low[choose]; s <= high[choose]; ++s) {
        if (dp[idx(choose, s)] == 1) {
            answer += s;
        }
    }
    return answer;
}

bool run_checkpoints() {
    if (solve(5, 3) != 330LL) {
        std::cerr << "Checkpoint failed for n=5,choose=3" << '\n';
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

    std::cout << solve(options.n, options.choose) << '\n';
    return 0;
}
