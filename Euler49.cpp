#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

bool is_prime(const int n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (int p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

std::string sorted_digits(const int n) {
    std::string s = std::to_string(n);
    std::sort(s.begin(), s.end());
    return s;
}

long long solve() {
    std::unordered_map<std::string, std::vector<int>> groups;

    for (int p = 1000; p < 10000; ++p) {
        if (is_prime(p)) {
            groups[sorted_digits(p)].push_back(p);
        }
    }

    for (auto& [key, values] : groups) {
        (void)key;
        if (values.size() < 3U) {
            continue;
        }

        std::sort(values.begin(), values.end());

        for (std::size_t i = 0; i < values.size(); ++i) {
            for (std::size_t j = i + 1; j < values.size(); ++j) {
                const int a = values[i];
                const int b = values[j];
                const int c = 2 * b - a;

                if (c > 9999) {
                    continue;
                }
                if (!std::binary_search(values.begin(), values.end(), c)) {
                    continue;
                }
                if (a == 1487 && b == 4817) {
                    continue;
                }

                return std::stoll(std::to_string(a) + std::to_string(b) + std::to_string(c));
            }
        }
    }

    return 0LL;
}

bool run_checkpoints() {
    if (!is_prime(1487) || !is_prime(4817) || !is_prime(8147)) {
        std::cerr << "Checkpoint failed for example primes" << '\n';
        return false;
    }
    if (sorted_digits(1487) != sorted_digits(4817) || sorted_digits(1487) != sorted_digits(8147)) {
        std::cerr << "Checkpoint failed for permutation grouping" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
