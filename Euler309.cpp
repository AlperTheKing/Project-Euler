#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int limit = 1000000;
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
    return options.limit >= 3;
}

i64 solve(const int limit) {
    const int lmt = limit - 1;
    std::vector<std::vector<std::pair<int, int>>> by_width(static_cast<std::size_t>(lmt + 1));
    std::set<std::pair<std::pair<int, int>, int>> uniq;

    const int root = static_cast<int>(std::sqrt(static_cast<double>(lmt))) + 2;
    for (int i = 0; i <= root; ++i) {
        for (int j = 1; j < i; ++j) {
            if (std::gcd(i, j) != 1) {
                continue;
            }
            int a = i * i - j * j;
            int b = 2 * i * j;
            int c = i * i + j * j;
            if (a > b) {
                std::swap(a, b);
            }
            for (int scale = 1; scale * c <= lmt; ++scale) {
                by_width[static_cast<std::size_t>(a * scale)].push_back({b * scale, c * scale});
                by_width[static_cast<std::size_t>(b * scale)].push_back({a * scale, c * scale});
            }
        }
    }

    for (int w = 1; w <= lmt; ++w) {
        auto& vec = by_width[static_cast<std::size_t>(w)];
        std::sort(vec.begin(), vec.end());
        for (const auto& p : vec) {
            for (const auto& q : vec) {
                if (p.first == q.first) {
                    break;
                }
                if ((static_cast<i64>(p.first) * static_cast<i64>(q.first)) %
                        static_cast<i64>(p.first + q.first) == 0) {
                    uniq.insert({{p.first, q.first}, w});
                }
            }
        }
    }

    return static_cast<i64>(uniq.size());
}

bool run_checkpoints() {
    if (solve(200) != 5LL) {
        std::cerr << "Checkpoint failed for x<y<200" << '\n';
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
