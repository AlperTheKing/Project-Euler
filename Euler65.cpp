#include <boost/multiprecision/cpp_int.hpp>

#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int term = 100;
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
        if (parse_int_after_prefix(arg, "--term=", options.term)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.term >= 1;
}

int coefficient_for_e(const int index) {
    if (index == 0) {
        return 2;
    }
    if (index % 3 == 2) {
        return 2 * ((index + 1) / 3);
    }
    return 1;
}

cpp_int convergent_numerator(const int term) {
    cpp_int num_prev = 1;
    cpp_int num = coefficient_for_e(0);
    cpp_int den_prev = 0;
    cpp_int den = 1;

    if (term == 1) {
        return num;
    }

    for (int i = 1; i < term; ++i) {
        const int a = coefficient_for_e(i);
        const cpp_int next_num = static_cast<cpp_int>(a) * num + num_prev;
        const cpp_int next_den = static_cast<cpp_int>(a) * den + den_prev;
        num_prev = num;
        den_prev = den;
        num = next_num;
        den = next_den;
    }

    return num;
}

int digit_sum(const cpp_int& x) {
    const std::string s = x.convert_to<std::string>();
    int sum = 0;
    for (char c : s) {
        sum += c - '0';
    }
    return sum;
}

int solve(const int term) {
    return digit_sum(convergent_numerator(term));
}

bool run_checkpoints() {
    if (solve(10) != 17) {
        std::cerr << "Checkpoint failed for term=10" << '\n';
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

    std::cout << solve(options.term) << '\n';
    return 0;
}
