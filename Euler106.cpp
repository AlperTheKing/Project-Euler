#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int n = 12;
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.n >= 2;
}

std::vector<int> masks_of_size(const int n, const int k) {
    std::vector<int> out;
    const int total = 1 << n;
    for (int mask = 0; mask < total; ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) == k) {
            out.push_back(mask);
        }
    }
    return out;
}

bool auto_ordered(const int a_mask, const int b_mask, const int n) {
    std::vector<int> a;
    std::vector<int> b;

    for (int i = 0; i < n; ++i) {
        if (((a_mask >> i) & 1) != 0) {
            a.push_back(i);
        }
        if (((b_mask >> i) & 1) != 0) {
            b.push_back(i);
        }
    }

    bool a_all_less = true;
    bool b_all_less = true;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] >= b[i]) {
            a_all_less = false;
        }
        if (b[i] >= a[i]) {
            b_all_less = false;
        }
    }

    return a_all_less || b_all_less;
}

int solve(const int n) {
    int count = 0;

    for (int k = 2; k <= n / 2; ++k) {
        const std::vector<int> subsets = masks_of_size(n, k);

        for (std::size_t i = 0; i < subsets.size(); ++i) {
            for (std::size_t j = i + 1; j < subsets.size(); ++j) {
                const int a = subsets[i];
                const int b = subsets[j];
                if ((a & b) != 0) {
                    continue;
                }
                if (auto_ordered(a, b, n)) {
                    continue;
                }
                ++count;
            }
        }
    }

    return count;
}

bool run_checkpoints() {
    if (solve(4) != 1) {
        std::cerr << "Checkpoint failed for n=4" << '\n';
        return false;
    }
    if (solve(7) != 70) {
        std::cerr << "Checkpoint failed for n=7" << '\n';
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
