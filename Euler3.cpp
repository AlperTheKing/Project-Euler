#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 n = 600851475143ULL;
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

u64 solve(u64 n) {
    u64 largest = 0ULL;

    while ((n & 1ULL) == 0ULL) {
        largest = 2ULL;
        n >>= 1ULL;
    }

    for (u64 d = 3ULL; d <= n / d; d += 2ULL) {
        while (n % d == 0ULL) {
            largest = d;
            n /= d;
        }
    }

    if (n > 1ULL) {
        largest = n;
    }

    return largest;
}

bool run_checkpoints() {
    if (solve(13195ULL) != 29ULL) {
        std::cerr << "Checkpoint failed for n=13195" << '\n';
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
