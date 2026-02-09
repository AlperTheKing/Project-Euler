#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

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
    try {
        value = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2;
}

u64 solve(const int limit) {
    const int n_max = limit - 1;
    std::vector<int> phi(static_cast<std::size_t>(n_max + 1), 0);
    std::vector<int> primes;
    std::vector<bool> is_composite(static_cast<std::size_t>(n_max + 1), false);

    phi[1] = 1;
    for (int n = 2; n <= n_max; ++n) {
        if (!is_composite[static_cast<std::size_t>(n)]) {
            primes.push_back(n);
            phi[static_cast<std::size_t>(n)] = n - 1;
        }
        for (const int p : primes) {
            const i64 m = static_cast<i64>(n) * p;
            if (m > n_max) {
                break;
            }
            is_composite[static_cast<std::size_t>(m)] = true;
            if (n % p == 0) {
                phi[static_cast<std::size_t>(m)] = phi[static_cast<std::size_t>(n)] * p;
                break;
            }
            phi[static_cast<std::size_t>(m)] = phi[static_cast<std::size_t>(n)] * (p - 1);
        }
    }

    u64 sum = 0;
    for (int n = 1; n <= n_max; ++n) {
        if (std::gcd(n, phi[static_cast<std::size_t>(n)]) != 1) {
            sum += static_cast<u64>(n);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(10) != 27ULL) {
        std::cerr << "Checkpoint failed: limit=10\n";
        return false;
    }
    if (solve(100) != 3186ULL) {
        std::cerr << "Checkpoint failed: limit=100\n";
        return false;
    }
    if (solve(1000) != 340576ULL) {
        std::cerr << "Checkpoint failed: limit=1000\n";
        return false;
    }
    if (solve(10000) != 34533819ULL) {
        std::cerr << "Checkpoint failed: limit=10000\n";
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
