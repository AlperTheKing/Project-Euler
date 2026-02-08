#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int m = 100;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--m=", options.m)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.m >= 1;
}

u64 solve(const int m) {
    std::vector<std::vector<int>> g(static_cast<std::size_t>(m + 1), std::vector<int>(static_cast<std::size_t>(m + 1), 0));
    for (int a = 1; a <= m; ++a) {
        for (int b = 1; b <= m; ++b) {
            g[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = std::gcd(a, b);
        }
    }

    const int max_i = 2 * m * m + 1;
    std::vector<unsigned char> is_square(static_cast<std::size_t>(max_i + 1), 0U);
    for (int x = 0; x * x <= max_i; ++x) {
        is_square[static_cast<std::size_t>(x * x)] = 1U;
    }

    u64 count = 0ULL;
    for (int a = 1; a <= m; ++a) {
        for (int b = 1; b <= m; ++b) {
            const int gab = g[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
            for (int c = 1; c <= m; ++c) {
                const int gbc = g[static_cast<std::size_t>(b)][static_cast<std::size_t>(c)];
                const int ac = a + c;
                for (int d = 1; d <= m; ++d) {
                    const int numerator =
                        ac * (b + d) -
                        gab -
                        gbc -
                        g[static_cast<std::size_t>(c)][static_cast<std::size_t>(d)] -
                        g[static_cast<std::size_t>(d)][static_cast<std::size_t>(a)];
                    const int interior = numerator / 2 + 1;
                    if (is_square[static_cast<std::size_t>(interior)] != 0U) {
                        ++count;
                    }
                }
            }
        }
    }
    return count;
}

u64 brute_small(const int m) {
    u64 count = 0ULL;
    for (int a = 1; a <= m; ++a) {
        for (int b = 1; b <= m; ++b) {
            for (int c = 1; c <= m; ++c) {
                for (int d = 1; d <= m; ++d) {
                    const int area2 = a * b + b * c + c * d + d * a;
                    const int boundary = std::gcd(a, b) + std::gcd(b, c) + std::gcd(c, d) + std::gcd(d, a);
                    const int interior = area2 / 2 - boundary / 2 + 1;
                    const int root = static_cast<int>(std::sqrt(static_cast<double>(interior)));
                    if (root * root == interior) {
                        ++count;
                    }
                }
            }
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(4) != 42ULL) {
        std::cerr << "Checkpoint failed: m=4 sample should be 42" << '\n';
        return false;
    }
    if (solve(7) != brute_small(7)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for m=7" << '\n';
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
    std::cout << solve(options.m) << '\n';
    return 0;
}
