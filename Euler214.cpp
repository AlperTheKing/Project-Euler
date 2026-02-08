#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 40000000;
    int target_length = 25;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--target-length=", options.target_length)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2 && options.target_length >= 1;
}

u64 solve(const int limit, const int target_length) {
    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> phi(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / 10));

    phi[1] = 1;
    for (int i = 2; i <= limit; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            phi[static_cast<std::size_t>(i)] = i - 1;
        }

        for (int p : primes) {
            const std::int64_t v = static_cast<std::int64_t>(i) * p;
            if (v > limit) {
                break;
            }

            is_composite[static_cast<std::size_t>(v)] = 1;
            if ((i % p) == 0) {
                phi[static_cast<std::size_t>(v)] = phi[static_cast<std::size_t>(i)] * p;
                break;
            }
            phi[static_cast<std::size_t>(v)] = phi[static_cast<std::size_t>(i)] * (p - 1);
        }
    }

    std::vector<std::uint8_t> chain(static_cast<std::size_t>(limit + 1), 0);
    chain[1] = 1;
    for (int n = 2; n <= limit; ++n) {
        chain[static_cast<std::size_t>(n)] = static_cast<std::uint8_t>(chain[static_cast<std::size_t>(phi[static_cast<std::size_t>(n)])] + 1);
    }

    u64 sum = 0;
    for (int p : primes) {
        if (chain[static_cast<std::size_t>(p)] == target_length) {
            sum += static_cast<u64>(p);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(100, 4) != 12ULL) {
        std::cerr << "Checkpoint failed for limit=100,target=4" << '\n';
        return false;
    }
    if (solve(1000, 5) != 43ULL) {
        std::cerr << "Checkpoint failed for limit=1000,target=5" << '\n';
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

    std::cout << solve(options.limit, options.target_length) << '\n';
    return 0;
}
