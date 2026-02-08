#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 10000000;
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

bool is_permutation(int a, int b) {
    std::array<int, 10> cnt{};
    while (a > 0) {
        ++cnt[static_cast<std::size_t>(a % 10)];
        a /= 10;
    }
    while (b > 0) {
        --cnt[static_cast<std::size_t>(b % 10)];
        b /= 10;
    }
    for (int v : cnt) {
        if (v != 0) {
            return false;
        }
    }
    return true;
}

int solve(const int limit) {
    std::vector<int> phi(static_cast<std::size_t>(limit));
    for (int i = 0; i < limit; ++i) {
        phi[static_cast<std::size_t>(i)] = i;
    }

    for (int i = 2; i < limit; ++i) {
        if (phi[static_cast<std::size_t>(i)] == i) {
            for (int j = i; j < limit; j += i) {
                phi[static_cast<std::size_t>(j)] -= phi[static_cast<std::size_t>(j)] / i;
            }
        }
    }

    int best_n = -1;
    int best_phi = 1;

    for (int n = 2; n < limit; ++n) {
        const int ph = phi[static_cast<std::size_t>(n)];
        if (!is_permutation(n, ph)) {
            continue;
        }

        if (best_n == -1 || static_cast<long long>(n) * best_phi <
                                 static_cast<long long>(best_n) * ph) {
            best_n = n;
            best_phi = ph;
        }
    }

    return best_n;
}

bool run_checkpoints() {
    if (solve(100) != 21) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
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
