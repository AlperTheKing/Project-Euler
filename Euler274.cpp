#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int prime_limit = 10000000;
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
        if (parse_int_after_prefix(arg, "--prime-limit=", options.prime_limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.prime_limit >= 3;
}

std::vector<int> primes_up_to(const int n) {
    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;

    for (int i = 2; i <= n; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
        for (int p : primes) {
            const long long v = 1LL * i * p;
            if (v > n) {
                break;
            }
            is_composite[static_cast<std::size_t>(v)] = 1;
            if (i % p == 0) {
                break;
            }
        }
    }

    return primes;
}

int inverse_mod_10(const int p) {
    int a = 10;
    int b = p;
    int x0 = 1;
    int x1 = 0;

    while (b != 0) {
        const int q = a / b;
        const int na = b;
        const int nb = a - q * b;
        a = na;
        b = nb;

        const int nx = x1;
        const int ny = x0 - q * x1;
        x0 = nx;
        x1 = ny;
    }

    int inv = x0 % p;
    if (inv < 0) {
        inv += p;
    }
    return inv;
}

u64 solve(const int prime_limit) {
    const std::vector<int> primes = primes_up_to(prime_limit - 1);

    u64 sum = 0;
    for (int p : primes) {
        if (p == 2 || p == 5) {
            continue;
        }
        sum += static_cast<u64>(inverse_mod_10(p));
    }
    return sum;
}

bool run_checkpoints() {
    if (inverse_mod_10(113) != 34) {
        std::cerr << "Checkpoint failed for p=113" << '\n';
        return false;
    }
    if (solve(1000) != 39517ULL) {
        std::cerr << "Checkpoint failed for prime limit 1000" << '\n';
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

    std::cout << solve(options.prime_limit) << '\n';
    return 0;
}
