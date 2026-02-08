#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    u64 limit = (1ULL << 50);
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

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 1;
}

u64 solve(const u64 limit) {
    const int max_k = static_cast<int>(std::sqrt(static_cast<long double>(limit)));

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(max_k / 10));
    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(max_k + 1), 0);
    std::vector<std::int8_t> mu(static_cast<std::size_t>(max_k + 1), 0);

    mu[1] = 1;
    for (int i = 2; i <= max_k; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }

        for (int p : primes) {
            const i64 v = static_cast<i64>(i) * p;
            if (v > max_k) {
                break;
            }
            is_composite[static_cast<std::size_t>(v)] = 1;
            if ((i % p) == 0) {
                mu[static_cast<std::size_t>(v)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(v)] = static_cast<std::int8_t>(-mu[static_cast<std::size_t>(i)]);
        }
    }

    i64 count = 0;
    for (i64 k = 1; k <= max_k; ++k) {
        const std::int8_t m = mu[static_cast<std::size_t>(k)];
        if (m == 0) {
            continue;
        }
        count += static_cast<i64>(m) * static_cast<i64>(limit / static_cast<u64>(k * k));
    }

    return static_cast<u64>(count);
}

bool run_checkpoints() {
    if (solve(100ULL) != 61ULL) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
        return false;
    }
    if (solve(1024ULL) != 624ULL) {
        std::cerr << "Checkpoint failed for limit=1024" << '\n';
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
