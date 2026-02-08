#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int min_divisors_exclusive = 500;
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
        if (parse_int_after_prefix(arg, "--divisors=", options.min_divisors_exclusive)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.min_divisors_exclusive >= 0;
}

u64 divisor_count(u64 n) {
    if (n == 1ULL) {
        return 1ULL;
    }

    u64 count = 1ULL;

    int exponent = 0;
    while ((n & 1ULL) == 0ULL) {
        n >>= 1ULL;
        ++exponent;
    }
    if (exponent > 0) {
        count *= static_cast<u64>(exponent + 1);
    }

    for (u64 p = 3ULL; p <= n / p; p += 2ULL) {
        exponent = 0;
        while (n % p == 0ULL) {
            n /= p;
            ++exponent;
        }
        if (exponent > 0) {
            count *= static_cast<u64>(exponent + 1);
        }
    }

    if (n > 1ULL) {
        count *= 2ULL;
    }

    return count;
}

u64 solve(const int min_divisors_exclusive) {
    for (u64 n = 1ULL;; ++n) {
        u64 a = n;
        u64 b = n + 1ULL;
        if ((a & 1ULL) == 0ULL) {
            a >>= 1ULL;
        } else {
            b >>= 1ULL;
        }

        const u64 divisors = divisor_count(a) * divisor_count(b);
        if (divisors > static_cast<u64>(min_divisors_exclusive)) {
            return n * (n + 1ULL) / 2ULL;
        }
    }
}

bool run_checkpoints() {
    if (solve(5) != 28ULL) {
        std::cerr << "Checkpoint failed for divisors>5" << '\n';
        return false;
    }
    if (solve(1) != 3ULL) {
        std::cerr << "Checkpoint failed for divisors>1" << '\n';
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

    std::cout << solve(options.min_divisors_exclusive) << '\n';
    return 0;
}
