#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int max_sum = 2011;
    int target_nines = 2011;
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
        if (parse_int_after_prefix(arg, "--max-sum=", options.max_sum) ||
            parse_int_after_prefix(arg, "--target-nines=", options.target_nines)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_sum >= 3 && options.target_nines >= 1;
}

u64 minimal_n(const int p, const int q, const int target_nines) {
    const long double beta = std::sqrt(static_cast<long double>(q)) -
                             std::sqrt(static_cast<long double>(p));
    const long double beta2 = beta * beta;
    const long double lambda = -std::log10(beta2);
    const long double target = static_cast<long double>(target_nines);

    u64 n = static_cast<u64>(target / lambda);
    if (n == 0ULL) {
        n = 1ULL;
    }
    while (static_cast<long double>(n) * lambda < target - 1e-15L) {
        ++n;
    }
    while (n > 1ULL && static_cast<long double>(n - 1ULL) * lambda >= target - 1e-15L) {
        --n;
    }
    return n;
}

u64 solve(const int max_sum, const int target_nines) {
    u64 total = 0ULL;
    for (int p = 1; p < max_sum; ++p) {
        for (int q = p + 1; p + q <= max_sum; ++q) {
            const long double beta = std::sqrt(static_cast<long double>(q)) -
                                     std::sqrt(static_cast<long double>(p));
            if (beta >= 1.0L) {
                continue;
            }
            total += minimal_n(p, q, target_nines);
        }
    }
    return total;
}

bool run_checkpoints() {
    if (minimal_n(2, 3, 1) != 2ULL || minimal_n(2, 3, 2) != 3ULL || minimal_n(2, 3, 3) != 4ULL) {
        std::cerr << "Checkpoint failed for sqrt(2)+sqrt(3) nines progression" << '\n';
        return false;
    }
    if (solve(5, 1) != 8ULL) {
        std::cerr << "Checkpoint failed for small max-sum=5 target=1" << '\n';
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
    std::cout << solve(options.max_sum, options.target_nines) << '\n';
    return 0;
}
