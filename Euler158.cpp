#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int alphabet = 26;
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
        if (parse_int_after_prefix(arg, "--alphabet=", options.alphabet)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.alphabet >= 1;
}

u64 binom(const int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);

    u64 result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * static_cast<u64>(n - k + i) / static_cast<u64>(i);
    }
    return result;
}

u64 pow2(const int exp) {
    return 1ULL << exp;
}

u64 p_of_n(const int alphabet, const int n) {
    const u64 eulerian_n1 = pow2(n) - static_cast<u64>(n) - 1ULL;
    return binom(alphabet, n) * eulerian_n1;
}

u64 solve(const int alphabet) {
    u64 best = 0;
    for (int n = 1; n <= alphabet; ++n) {
        best = std::max(best, p_of_n(alphabet, n));
    }
    return best;
}

bool run_checkpoints() {
    if (p_of_n(26, 3) != 10400ULL) {
        std::cerr << "Checkpoint failed for p(3)" << '\n';
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

    std::cout << solve(options.alphabet) << '\n';
    return 0;
}
