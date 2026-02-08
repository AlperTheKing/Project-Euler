#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int limit = 100000;
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
    return options.limit >= 0;
}

std::vector<int> build_grundy(const int limit) {
    std::vector<int> g(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> seen(512, -1);

    for (int n = 1; n <= limit; ++n) {
        for (int k = 1; k * k <= n; ++k) {
            const int val = g[static_cast<std::size_t>(n - k * k)];
            if (val >= static_cast<int>(seen.size())) {
                seen.resize(static_cast<std::size_t>(val + 64), -1);
            }
            seen[static_cast<std::size_t>(val)] = n;
        }
        while (seen[static_cast<std::size_t>(g[static_cast<std::size_t>(n)])] == n) {
            ++g[static_cast<std::size_t>(n)];
        }
    }

    return g;
}

i64 solve(const int limit) {
    const std::vector<int> g = build_grundy(limit);
    int max_g = 0;
    for (int x : g) {
        max_g = std::max(max_g, x);
    }

    std::vector<i64> prefix(static_cast<std::size_t>(max_g + 1), 0);
    std::vector<i64> suffix(static_cast<std::size_t>(max_g + 1), 0);
    for (int n = 0; n <= limit; ++n) {
        ++suffix[static_cast<std::size_t>(g[static_cast<std::size_t>(n)])];
    }

    i64 answer = 0;
    for (int b = 0; b <= limit; ++b) {
        const int gb = g[static_cast<std::size_t>(b)];
        ++prefix[static_cast<std::size_t>(gb)];

        for (int ga = 0; ga <= max_g; ++ga) {
            const int gc = ga ^ gb;
            if (gc > max_g) {
                continue;
            }
            answer += prefix[static_cast<std::size_t>(ga)] * suffix[static_cast<std::size_t>(gc)];
        }

        --suffix[static_cast<std::size_t>(gb)];
    }

    return answer;
}

i64 brute_small(const int limit) {
    const std::vector<int> g = build_grundy(limit);
    i64 count = 0;
    for (int a = 0; a <= limit; ++a) {
        for (int b = a; b <= limit; ++b) {
            for (int c = b; c <= limit; ++c) {
                if ((g[static_cast<std::size_t>(a)] ^ g[static_cast<std::size_t>(b)] ^ g[static_cast<std::size_t>(c)]) == 0) {
                    ++count;
                }
            }
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(29) != 1160LL) {
        std::cerr << "Checkpoint failed for limit 29" << '\n';
        return false;
    }
    if (solve(40) != brute_small(40)) {
        std::cerr << "Checkpoint failed for brute cross-check at limit 40" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
