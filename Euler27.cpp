#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int abs_a = 999;
    int abs_b = 999;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--abs-a=", options.abs_a)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--abs-b=", options.abs_b)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.abs_a >= 0 && options.abs_b >= 0;
}

bool is_prime(const i64 n) {
    if (n < 2) {
        return false;
    }
    if ((n & 1LL) == 0LL) {
        return n == 2;
    }
    for (i64 p = 3; p <= n / p; p += 2) {
        if (n % p == 0LL) {
            return false;
        }
    }
    return true;
}

int consecutive_prime_count(const int a, const int b) {
    int n = 0;
    while (true) {
        const i64 value = static_cast<i64>(n) * static_cast<i64>(n) +
                          static_cast<i64>(a) * static_cast<i64>(n) +
                          static_cast<i64>(b);
        if (!is_prime(value)) {
            break;
        }
        ++n;
    }
    return n;
}

int solve(const int abs_a, const int abs_b) {
    int best_a = 0;
    int best_b = 0;
    int best_count = -1;

    for (int a = -abs_a; a <= abs_a; ++a) {
        for (int b = -abs_b; b <= abs_b; ++b) {
            if (!is_prime(b)) {
                continue;
            }
            const int count = consecutive_prime_count(a, b);
            if (count > best_count) {
                best_count = count;
                best_a = a;
                best_b = b;
            }
        }
    }

    return best_a * best_b;
}

bool run_checkpoints() {
    if (consecutive_prime_count(1, 41) != 40) {
        std::cerr << "Checkpoint failed for n^2+n+41" << '\n';
        return false;
    }
    if (consecutive_prime_count(-79, 1601) != 80) {
        std::cerr << "Checkpoint failed for n^2-79n+1601" << '\n';
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

    std::cout << solve(options.abs_a, options.abs_b) << '\n';
    return 0;
}
