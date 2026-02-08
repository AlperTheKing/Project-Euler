#include <boost/multiprecision/cpp_int.hpp>

#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int expansions = 1000;
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
        if (parse_int_after_prefix(arg, "--expansions=", options.expansions)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.expansions >= 1;
}

int digit_count(const cpp_int& x) {
    return static_cast<int>(x.convert_to<std::string>().size());
}

int solve(const int expansions) {
    cpp_int num = 3;  // first expansion numerator
    cpp_int den = 2;  // first expansion denominator

    int count = 0;

    for (int i = 1; i <= expansions; ++i) {
        if (digit_count(num) > digit_count(den)) {
            ++count;
        }

        const cpp_int next_num = num + 2 * den;
        const cpp_int next_den = num + den;
        num = next_num;
        den = next_den;
    }

    return count;
}

bool run_checkpoints() {
    if (solve(8) != 1) {
        std::cerr << "Checkpoint failed for first 8 expansions" << '\n';
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

    std::cout << solve(options.expansions) << '\n';
    return 0;
}
