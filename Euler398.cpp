#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    long long n = 10000000;
    int m = 100;
    bool run_checkpoints = true;
};

bool parse_ll_after_prefix(const std::string& arg, const std::string& prefix, long long& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoll(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_ll_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--m=", options.m)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.m >= 2 && options.m <= options.n;
}

long double ratio_binom_drop(long long a, int drop, int choose_r) {
    // Computes C(a-drop, r) / C(a, r) by telescoping product.
    long double ratio = 1.0L;
    for (int j = 0; j < choose_r; ++j) {
        ratio *= static_cast<long double>(a - drop - j) / static_cast<long double>(a - j);
    }
    return ratio;
}

long double expected_second_shortest(long long n, int m) {
    // P(X_(2) >= k) = [ m*C(n-(m-1)(k-1)-1,m-1) - (m-1)*C(n-m(k-1)-1,m-1) ] / C(n-1,m-1).
    // Then E[X_(2)] = sum_{k>=1} P(X_(2) >= k).
    const int r = m - 1;

    long long a_r = n - 1;  // numerator top for first binomial ratio.
    long long a_s = n - 1;  // numerator top for second binomial ratio.
    long double ratio_r = 1.0L;
    long double ratio_s = 1.0L;

    long double answer = 0.0L;
    while (true) {
        const long double tail = static_cast<long double>(m) * ratio_r -
                                 static_cast<long double>(m - 1) * ratio_s;
        answer += tail;

        bool has_next = false;
        if (a_r - (m - 1) >= r) {
            ratio_r *= ratio_binom_drop(a_r, m - 1, r);
            a_r -= (m - 1);
            has_next = true;
        } else {
            ratio_r = 0.0L;
        }

        if (a_s - m >= r) {
            ratio_s *= ratio_binom_drop(a_s, m, r);
            a_s -= m;
            has_next = true;
        } else {
            ratio_s = 0.0L;
        }

        if (!has_next) {
            break;
        }
    }

    return answer;
}

bool nearly_equal(long double a, long double b, long double rel_tol = 1e-12L) {
    const long double scale = std::max(std::fabsl(a), std::fabsl(b));
    if (scale == 0.0L) {
        return true;
    }
    return std::fabsl(a - b) <= rel_tol * scale;
}

bool run_checkpoints() {
    if (!nearly_equal(expected_second_shortest(3, 2), 2.0L, 1e-14L)) {
        std::cerr << "Checkpoint failed: E(3,2)\n";
        return false;
    }

    const long double e83 = expected_second_shortest(8, 3);
    if (!nearly_equal(e83, 16.0L / 7.0L, 1e-13L)) {
        std::cerr << "Checkpoint failed: E(8,3)\n";
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

    const long double answer = expected_second_shortest(options.n, options.m);
    std::cout << std::fixed << std::setprecision(5) << static_cast<double>(answer) << '\n';
    return 0;
}
