#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int from = 5;
    int to = 10000;
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
        if (parse_int_after_prefix(arg, "--from=", options.from) ||
            parse_int_after_prefix(arg, "--to=", options.to)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.from >= 2 && options.to >= options.from;
}

int best_k(const int n) {
    int k = static_cast<int>(std::floor(static_cast<double>(n) / std::exp(1.0)));
    if (k < 1) {
        k = 1;
    }

    auto score = [&](const int x) {
        return static_cast<double>(x) * (std::log(static_cast<double>(n)) - std::log(static_cast<double>(x)));
    };

    const int k2 = k + 1;
    return (score(k2) > score(k)) ? k2 : k;
}

int D(const int n) {
    const int k = best_k(n);
    int q = k / std::gcd(n, k);

    while ((q % 2) == 0) {
        q /= 2;
    }
    while ((q % 5) == 0) {
        q /= 5;
    }

    return (q == 1) ? -n : n;
}

i64 solve(const int from, const int to) {
    i64 sum = 0;
    for (int n = from; n <= to; ++n) {
        sum += D(n);
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(5, 100) != 2438) {
        std::cerr << "Checkpoint failed for [5,100]" << '\n';
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

    std::cout << solve(options.from, options.to) << '\n';
    return 0;
}
