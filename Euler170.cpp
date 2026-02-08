#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

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

bool is_pandigital_0_to_9(const std::string& s) {
    if (s.size() != 10U) {
        return false;
    }
    std::array<int, 10> freq{};
    for (char c : s) {
        if (c < '0' || c > '9') {
            return false;
        }
        ++freq[static_cast<std::size_t>(c - '0')];
    }
    for (int f : freq) {
        if (f != 1) {
            return false;
        }
    }
    return true;
}

bool is_pandigital_1_to_9(const std::string& s) {
    if (s.size() != 9U) {
        return false;
    }
    std::array<int, 10> freq{};
    for (char c : s) {
        if (c < '1' || c > '9') {
            return false;
        }
        ++freq[static_cast<std::size_t>(c - '0')];
    }
    for (int d = 1; d <= 9; ++d) {
        if (freq[static_cast<std::size_t>(d)] != 1) {
            return false;
        }
    }
    return true;
}

void collect_divisors(const int n, std::vector<int>& divisors) {
    divisors.clear();
    for (int d = 1; static_cast<long long>(d) * d <= n; ++d) {
        if (n % d != 0) {
            continue;
        }
        divisors.push_back(d);
        if (d * d != n) {
            divisors.push_back(n / d);
        }
    }
    std::sort(divisors.begin(), divisors.end(), std::greater<int>());
}

bool output_has_valid_input(const std::string& output) {
    std::vector<int> divisors;

    for (int split_mask = 1; split_mask < (1 << 9); ++split_mask) {
        std::vector<int> products;
        int start = 0;
        bool good_partition = true;
        for (int bit = 0; bit < 9; ++bit) {
            if (((split_mask >> bit) & 1) == 0) {
                continue;
            }
            const int end = bit + 1;
            if (output[static_cast<std::size_t>(start)] == '0' && (end - start) > 1) {
                good_partition = false;
                break;
            }
            products.push_back(std::stoi(output.substr(static_cast<std::size_t>(start),
                                                       static_cast<std::size_t>(end - start))));
            start = end;
        }
        if (!good_partition) {
            continue;
        }
        if (output[static_cast<std::size_t>(start)] == '0' &&
            (static_cast<int>(output.size()) - start) > 1) {
            continue;
        }
        products.push_back(std::stoi(output.substr(static_cast<std::size_t>(start))));
        if (products.size() < 2U) {
            continue;
        }

        int g = products.front();
        for (int value : products) {
            g = std::gcd(g, value);
        }

        collect_divisors(g, divisors);
        for (int multiplier : divisors) {
            if (multiplier == 0) {
                continue;
            }

            std::string input_concat = std::to_string(multiplier);
            bool divisible = true;
            for (int value : products) {
                if (value % multiplier != 0) {
                    divisible = false;
                    break;
                }
                input_concat += std::to_string(value / multiplier);
            }
            if (!divisible) {
                continue;
            }
            if (is_pandigital_0_to_9(input_concat)) {
                return true;
            }
        }
    }

    return false;
}

u64 solve() {
    std::array<int, 10> digits{9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    do {
        if (digits[0] == 0) {
            continue;
        }
        std::string output;
        output.reserve(10);
        for (int d : digits) {
            output.push_back(static_cast<char>('0' + d));
        }

        if (output_has_valid_input(output)) {
            return static_cast<u64>(std::stoull(output));
        }
    } while (std::prev_permutation(digits.begin(), digits.end()));

    return 0;
}

bool run_checkpoints() {
    const int k = 6;
    const int a = 1273;
    const int b = 9854;
    const std::string product_concat = std::to_string(k * a) + std::to_string(k * b);
    const std::string input_concat = std::to_string(k) + std::to_string(a) + std::to_string(b);
    if (!is_pandigital_1_to_9(product_concat) || !is_pandigital_1_to_9(input_concat)) {
        std::cerr << "Checkpoint failed for 1-to-9 sample" << '\n';
        return false;
    }

    const u64 ans = solve();
    const std::string ans_s = std::to_string(ans);
    if (!is_pandigital_0_to_9(ans_s) || !output_has_valid_input(ans_s)) {
        std::cerr << "Checkpoint failed for final candidate validation" << '\n';
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
