#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 n = 1'000'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value_out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value_out = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
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
    return options.n >= 4ULL;
}

u64 smallest_prime_factor(const u64 n) {
    if ((n & 1ULL) == 0ULL) {
        return 2ULL;
    }
    if (n % 3ULL == 0ULL) {
        return 3ULL;
    }
    for (u64 d = 5ULL; d <= n / d; d += 6ULL) {
        if (n % d == 0ULL) {
            return d;
        }
        if (n % (d + 2ULL) == 0ULL) {
            return d + 2ULL;
        }
    }
    return n;
}

u64 brute_g(const u64 n) {
    u64 g = 13ULL;
    for (u64 k = 5ULL; k <= n; ++k) {
        g += std::gcd(k, g);
    }
    return g;
}

u64 solve(const u64 n) {
    if (n <= 9ULL) {
        return brute_g(n);
    }

    // First index where g(n)=3n is n=9.
    u64 anchor = 9ULL;

    while (anchor < n) {
        const u64 p = smallest_prime_factor(2ULL * anchor - 1ULL);
        const u64 next_anchor = anchor + (p - 1ULL) / 2ULL;
        if (next_anchor > n) {
            return n + 2ULL * anchor;
        }
        anchor = next_anchor;
    }

    return 3ULL * n;
}

bool run_checkpoints() {
    for (u64 n = 4ULL; n <= 5'000ULL; ++n) {
        if (solve(n) != brute_g(n)) {
            std::cerr << "Checkpoint failed: brute-force mismatch at n=" << n << '\n';
            return false;
        }
    }

    if (solve(1'000ULL) != 2'524ULL) {
        std::cerr << "Checkpoint failed: g(1000)\n";
        return false;
    }
    if (solve(1'000'000ULL) != 2'624'152ULL) {
        std::cerr << "Checkpoint failed: g(1000000)\n";
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
