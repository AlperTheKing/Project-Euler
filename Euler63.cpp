#include <boost/multiprecision/cpp_int.hpp>

#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_int;

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

int digit_count(const cpp_int& x) {
    return static_cast<int>(x.convert_to<std::string>().size());
}

int solve() {
    int total = 0;

    for (int base = 1; base <= 9; ++base) {
        cpp_int value = 1;
        for (int power = 1;; ++power) {
            value *= base;
            const int digits = digit_count(value);
            if (digits == power) {
                ++total;
            }
            if (digits < power) {
                break;
            }
        }
    }

    return total;
}

bool run_checkpoints() {
    cpp_int value = 1;
    for (int i = 0; i < 21; ++i) {
        value *= 9;
    }
    if (digit_count(value) != 21) {
        std::cerr << "Checkpoint failed for 9^21" << '\n';
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
