#include <boost/multiprecision/cpp_int.hpp>

#include <cmath>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int limit = 1000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2;
}

cpp_int minimal_x_for_pell(const int D) {
    const int a0 = static_cast<int>(std::sqrt(D));
    if (a0 * a0 == D) {
        return 0;
    }

    int m = 0;
    int d = 1;
    int a = a0;

    cpp_int num_prev = 1;
    cpp_int num = a;
    cpp_int den_prev = 0;
    cpp_int den = 1;

    while (num * num - static_cast<cpp_int>(D) * den * den != 1) {
        m = d * a - m;
        d = (D - m * m) / d;
        a = (a0 + m) / d;

        const cpp_int next_num = static_cast<cpp_int>(a) * num + num_prev;
        const cpp_int next_den = static_cast<cpp_int>(a) * den + den_prev;
        num_prev = num;
        den_prev = den;
        num = next_num;
        den = next_den;
    }

    return num;
}

int solve(const int limit) {
    int best_D = 0;
    cpp_int best_x = 0;

    for (int D = 2; D <= limit; ++D) {
        const cpp_int x = minimal_x_for_pell(D);
        if (x > best_x) {
            best_x = x;
            best_D = D;
        }
    }

    return best_D;
}

bool run_checkpoints() {
    if (solve(7) != 5) {
        std::cerr << "Checkpoint failed for limit=7" << '\n';
        return false;
    }
    if (solve(13) != 13) {
        std::cerr << "Checkpoint failed for limit=13" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
