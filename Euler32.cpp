#include <array>
#include <iostream>
#include <set>
#include <string>

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

bool is_1_to_9_pandigital_triplet(const int a, const int b, const int product) {
    std::array<bool, 10> used{};
    used.fill(false);

    int digits_seen = 0;
    auto add_digits = [&](int value) -> bool {
        while (value > 0) {
            const int d = value % 10;
            if (d == 0 || used[static_cast<std::size_t>(d)]) {
                return false;
            }
            used[static_cast<std::size_t>(d)] = true;
            ++digits_seen;
            value /= 10;
        }
        return true;
    };

    if (!add_digits(a) || !add_digits(b) || !add_digits(product)) {
        return false;
    }

    if (digits_seen != 9) {
        return false;
    }
    for (int d = 1; d <= 9; ++d) {
        if (!used[static_cast<std::size_t>(d)]) {
            return false;
        }
    }
    return true;
}

int solve() {
    std::set<int> products;

    for (int a = 1; a <= 9; ++a) {
        for (int b = 1234; b <= 9876; ++b) {
            const int p = a * b;
            if (is_1_to_9_pandigital_triplet(a, b, p)) {
                products.insert(p);
            }
        }
    }

    for (int a = 12; a <= 98; ++a) {
        for (int b = 123; b <= 987; ++b) {
            const int p = a * b;
            if (is_1_to_9_pandigital_triplet(a, b, p)) {
                products.insert(p);
            }
        }
    }

    int total = 0;
    for (const int p : products) {
        total += p;
    }
    return total;
}

bool run_checkpoints() {
    if (!is_1_to_9_pandigital_triplet(39, 186, 7254)) {
        std::cerr << "Checkpoint failed for example identity" << '\n';
        return false;
    }
    if (is_1_to_9_pandigital_triplet(12, 34, 408)) {
        std::cerr << "Checkpoint failed for non-pandigital identity" << '\n';
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
