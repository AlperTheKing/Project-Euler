#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int count = 25;
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
        if (parse_int_after_prefix(arg, "--count=", options.count)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.count >= 1;
}

bool is_prime(const int n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (int p = 3; static_cast<i64>(p) * p <= n; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

int A_of_n(const int n) {
    int rem = 1 % n;
    int k = 1;
    while (rem != 0) {
        rem = (rem * 10 + 1) % n;
        ++k;
    }
    return k;
}

i64 solve(const int wanted_count) {
    int found = 0;
    i64 sum = 0;

    for (int n = 3; found < wanted_count; n += 2) {
        if ((n % 5) == 0) {
            continue;
        }
        if (is_prime(n)) {
            continue;
        }
        if (std::gcd(n, 10) != 1) {
            continue;
        }

        const int a = A_of_n(n);
        if (((n - 1) % a) == 0) {
            ++found;
            sum += n;
        }
    }

    return sum;
}

bool run_checkpoints() {
    if (solve(5) != 1985) {
        std::cerr << "Checkpoint failed for first five values" << '\n';
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

    std::cout << solve(options.count) << '\n';
    return 0;
}
