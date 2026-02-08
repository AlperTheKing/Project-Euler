#include <boost/multiprecision/cpp_int.hpp>

#include <iostream>
#include <set>
#include <string>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int max_a = 100;
    int max_b = 100;
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
    for (const char c : tail) {
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
        if (parse_int_after_prefix(arg, "--max-a=", options.max_a)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--max-b=", options.max_b)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_a >= 2 && options.max_b >= 2;
}

std::size_t solve(const int max_a, const int max_b) {
    std::set<cpp_int> values;

    for (int a = 2; a <= max_a; ++a) {
        cpp_int value = 1;
        for (int b = 1; b <= max_b; ++b) {
            value *= a;
            if (b >= 2) {
                values.insert(value);
            }
        }
    }

    return values.size();
}

bool run_checkpoints() {
    if (solve(5, 5) != 15U) {
        std::cerr << "Checkpoint failed for range 2..5" << '\n';
        return false;
    }
    if (solve(10, 10) != 69U) {
        std::cerr << "Checkpoint failed for range 2..10" << '\n';
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

    std::cout << solve(options.max_a, options.max_b) << '\n';
    return 0;
}
