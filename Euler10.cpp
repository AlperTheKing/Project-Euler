#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 2000000;
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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
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

u64 solve(const int limit) {
    if (limit <= 2) {
        return 0ULL;
    }

    std::vector<bool> is_prime(static_cast<std::size_t>(limit), true);
    is_prime[0] = false;
    is_prime[1] = false;

    for (int p = 2; p * p < limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q < limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    u64 sum = 0ULL;
    for (int p = 2; p < limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            sum += static_cast<u64>(p);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(10) != 17ULL) {
        std::cerr << "Checkpoint failed for limit=10" << '\n';
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
