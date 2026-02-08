#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int family_size = 5;
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
        if (parse_int_after_prefix(arg, "--family-size=", options.family_size)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.family_size >= 2;
}

i64 solve(const int family_size) {
    std::unordered_map<std::string, std::vector<i64>> groups;
    int current_digits = 1;

    for (i64 n = 1;; ++n) {
        const i64 cube = n * n * n;
        const int digits = static_cast<int>(std::to_string(cube).size());

        if (digits != current_digits) {
            i64 best = -1;
            for (const auto& [key, vec] : groups) {
                (void)key;
                if (static_cast<int>(vec.size()) == family_size) {
                    const i64 candidate = *std::min_element(vec.begin(), vec.end());
                    if (best == -1 || candidate < best) {
                        best = candidate;
                    }
                }
            }
            if (best != -1) {
                return best;
            }
            groups.clear();
            current_digits = digits;
        }

        std::string key = std::to_string(cube);
        std::sort(key.begin(), key.end());
        groups[key].push_back(cube);
    }
}

bool run_checkpoints() {
    if (solve(3) != 41063625LL) {
        std::cerr << "Checkpoint failed for family size 3" << '\n';
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

    std::cout << solve(options.family_size) << '\n';
    return 0;
}
