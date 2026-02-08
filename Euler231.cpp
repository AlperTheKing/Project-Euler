#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int n = 20000000;
    int k = 15000000;
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
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--k=", options.k)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.k >= 0 && options.k <= options.n;
}

std::vector<int> primes_up_to(const int n) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(n + 1), 1);
    is_prime[0] = 0;
    is_prime[1] = 0;

    for (int i = 2; i * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

u64 exponent_in_factorial(const int n, const int p) {
    u64 exp = 0;
    int x = n;
    while (x > 0) {
        x /= p;
        exp += static_cast<u64>(x);
    }
    return exp;
}

u64 solve(const int n, const int k) {
    const int r = n - k;
    const std::vector<int> primes = primes_up_to(n);

    u64 sum = 0;
    for (int p : primes) {
        const u64 exp = exponent_in_factorial(n, p) -
                        exponent_in_factorial(k, p) -
                        exponent_in_factorial(r, p);
        sum += static_cast<u64>(p) * exp;
    }

    return sum;
}

bool run_checkpoints() {
    if (solve(10, 3) != 14ULL) {
        std::cerr << "Checkpoint failed for C(10,3)" << '\n';
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

    std::cout << solve(options.n, options.k) << '\n';
    return 0;
}
