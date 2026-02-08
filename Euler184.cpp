#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int radius = 105;
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
        if (parse_int_after_prefix(arg, "--radius=", options.radius)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.radius >= 2;
}

i64 solve(const int radius) {
    constexpr double kPi = 3.141592653589793238462643383279502884;
    constexpr double kEps = 1e-8;

    std::vector<double> angles;
    for (int x = -radius; x <= radius; ++x) {
        for (int y = -radius; y <= radius; ++y) {
            if ((x == 0 && y == 0) || x * x + y * y >= radius * radius) {
                continue;
            }
            angles.push_back(std::atan2(static_cast<double>(y), static_cast<double>(x)));
        }
    }

    std::sort(angles.begin(), angles.end());
    const int total = static_cast<int>(angles.size());
    if (total < 3) {
        return 0;
    }

    angles.reserve(static_cast<std::size_t>(2 * total));
    for (int i = 0; i < total; ++i) {
        angles.push_back(angles[static_cast<std::size_t>(i)] + 2.0 * kPi);
    }

    const i64 n = total;
    const i64 half = n / 2;
    i64 ans = n * (n - 1) * (n - 2) / 6;

    for (i64 i = 0; i < half;) {
        i64 cnt = 1;
        while (i + cnt < half &&
               angles[static_cast<std::size_t>(i + cnt)] - angles[static_cast<std::size_t>(i)] <= kEps) {
            ++cnt;
        }

        ans -= (half - cnt) * (half - cnt - 1) * cnt;
        ans -= cnt * cnt * (n - 2 * cnt) + cnt * (cnt - 1) * (half - cnt);

        const i64 doubled = 2 * cnt;
        ans -= doubled * (doubled - 1) * (doubled - 2) / 6;
        i += cnt;
    }

    return ans;
}

bool run_checkpoints() {
    if (solve(2) != 8LL) {
        std::cerr << "Checkpoint failed for r=2" << '\n';
        return false;
    }
    if (solve(3) != 360LL) {
        std::cerr << "Checkpoint failed for r=3" << '\n';
        return false;
    }
    if (solve(5) != 10600LL) {
        std::cerr << "Checkpoint failed for r=5" << '\n';
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

    std::cout << solve(options.radius) << '\n';
    return 0;
}
