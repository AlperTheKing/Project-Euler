#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int from_n = 1864;
    int to_n = 1909;
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
        if (parse_int_after_prefix(arg, "--from=", options.from_n) ||
            parse_int_after_prefix(arg, "--to=", options.to_n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.from_n >= 3 && options.to_n >= options.from_n;
}

std::vector<int> compute_totients(const int limit) {
    std::vector<int> phi(static_cast<std::size_t>(limit + 1));
    for (int i = 0; i <= limit; ++i) {
        phi[static_cast<std::size_t>(i)] = i;
    }
    for (int p = 2; p <= limit; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (int k = p; k <= limit; k += p) {
            phi[static_cast<std::size_t>(k)] -= phi[static_cast<std::size_t>(k)] / p;
        }
    }
    return phi;
}

u64 solve(const int from_n, const int to_n) {
    const std::vector<int> phi = compute_totients(to_n);
    std::vector<std::uint8_t> denominator_used(static_cast<std::size_t>(to_n + 1), 0U);

    for (int n = from_n; n <= to_n; ++n) {
        for (int d = 1; static_cast<long long>(d) * d <= n; ++d) {
            if (n % d != 0) {
                continue;
            }
            denominator_used[static_cast<std::size_t>(d)] = 1U;
            const int e = n / d;
            denominator_used[static_cast<std::size_t>(e)] = 1U;
        }
    }

    u64 sides = 0;
    for (int d = 1; d <= to_n; ++d) {
        if (denominator_used[static_cast<std::size_t>(d)] != 0U) {
            sides += static_cast<u64>(phi[static_cast<std::size_t>(d)]);
        }
    }
    return sides;
}

bool run_checkpoints() {
    if (solve(3, 4) != 6ULL) {
        std::cerr << "Checkpoint failed for S3 + S4 sample" << '\n';
        return false;
    }
    if (solve(4, 4) != 4ULL) {
        std::cerr << "Checkpoint failed for single polygon S4" << '\n';
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
    std::cout << solve(options.from_n, options.to_n) << '\n';
    return 0;
}
