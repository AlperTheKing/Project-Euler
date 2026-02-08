#include <iostream>
#include <numeric>
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

bool is_nontrivial_digit_canceling(const int num, const int den) {
    if (num >= den || num < 10 || den < 10) {
        return false;
    }

    const int n1 = num / 10;
    const int n2 = num % 10;
    const int d1 = den / 10;
    const int d2 = den % 10;

    if (n2 == 0 && d2 == 0) {
        return false;
    }

    auto equivalent = [&](int a, int b) -> bool {
        if (b == 0) {
            return false;
        }
        return static_cast<long long>(num) * b == static_cast<long long>(den) * a;
    };

    if (n1 == d1 && n1 != 0 && equivalent(n2, d2)) {
        return true;
    }
    if (n1 == d2 && n1 != 0 && equivalent(n2, d1)) {
        return true;
    }
    if (n2 == d1 && n2 != 0 && equivalent(n1, d2)) {
        return true;
    }
    if (n2 == d2 && n2 != 0 && equivalent(n1, d1)) {
        return true;
    }

    return false;
}

int solve() {
    long long product_num = 1;
    long long product_den = 1;

    for (int num = 10; num < 100; ++num) {
        for (int den = num + 1; den < 100; ++den) {
            if (is_nontrivial_digit_canceling(num, den)) {
                product_num *= num;
                product_den *= den;
            }
        }
    }

    const long long g = std::gcd(product_num, product_den);
    return static_cast<int>(product_den / g);
}

bool run_checkpoints() {
    if (!is_nontrivial_digit_canceling(49, 98)) {
        std::cerr << "Checkpoint failed for 49/98" << '\n';
        return false;
    }
    if (is_nontrivial_digit_canceling(30, 50)) {
        std::cerr << "Checkpoint failed for trivial cancellation" << '\n';
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
