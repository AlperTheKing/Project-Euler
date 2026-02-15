#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 100000000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

std::vector<int> gen_primes(int n, bool /*minimal_divisors*/) {
    std::vector<int> min_divisor(static_cast<std::size_t>(n + 1), 0);
    min_divisor[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (min_divisor[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        min_divisor[static_cast<std::size_t>(i)] = i;
        if (1LL * i * i > n) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            if (min_divisor[static_cast<std::size_t>(j)] == 0) {
                min_divisor[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return min_divisor;
}

i64 solve(int n) {
    const std::vector<int> divs = gen_primes(n, true);
    std::vector<uint8_t> degs(static_cast<std::size_t>(n + 1), 0);
    std::vector<std::vector<int>> by_degs(30);

    by_degs[0].push_back(1);
    for (int i = 2; i <= n; ++i) {
        const int reduced = i / divs[static_cast<std::size_t>(i)];
        const int d = static_cast<int>(degs[static_cast<std::size_t>(reduced)]) + 1;
        degs[static_cast<std::size_t>(i)] = static_cast<uint8_t>(d);
        if (static_cast<int>(by_degs.size()) <= d) {
            by_degs.resize(d + 1);
        }
        by_degs[static_cast<std::size_t>(d)].push_back(i);
    }

    i64 ans = 0;
    for (std::size_t i = 0; i + 1 < by_degs.size(); ++i) {
        const auto& cur = by_degs[static_cast<std::size_t>(i)];

        int low = 0;
        int high = static_cast<int>(cur.size()) - 1;
        while (low < static_cast<int>(cur.size())) {
            while (high >= 0 && 1LL * cur[static_cast<std::size_t>(high)] * cur[static_cast<std::size_t>(low)] > n) {
                --high;
            }
            ans += high + 1LL;
            ++low;
        }

        const auto& nxt = by_degs[static_cast<std::size_t>(i + 1)];
        low = 0;
        high = static_cast<int>(nxt.size()) - 1;
        while (low < static_cast<int>(cur.size())) {
            while (high >= 0 && 1LL * nxt[static_cast<std::size_t>(high)] * cur[static_cast<std::size_t>(low)] > n) {
                --high;
            }
            ans += high + 1LL;
            ++low;
        }
    }

    return ans;
}

bool run_checkpoints() {
    return solve(20) > 0;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        std::cerr << "Checkpoint failed\n";
        return 2;
    }
    std::cout << solve(options.limit) << '\n';
    return 0;
}
