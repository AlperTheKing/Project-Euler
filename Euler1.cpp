#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 limit = 1000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

u64 sum_of_multiples_below(const u64 limit, const u64 divisor) {
    const u64 n = (limit - 1ULL) / divisor;
    return divisor * n * (n + 1ULL) / 2ULL;
}

u64 solve(const u64 limit) {
    if (limit == 0ULL) {
        return 0ULL;
    }
    return sum_of_multiples_below(limit, 3ULL) + sum_of_multiples_below(limit, 5ULL) -
           sum_of_multiples_below(limit, 15ULL);
}

bool run_checkpoints() {
    if (solve(10ULL) != 23ULL) {
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
