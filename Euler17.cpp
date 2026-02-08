#include <cstdint>
#include <iostream>
#include <string>

namespace {

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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
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

    return options.limit >= 1 && options.limit <= 1000;
}

int letters_for_number(int n) {
    static const int ones[] = {
        0, 3, 3, 5, 4, 4, 3, 5, 5, 4,
    };
    static const int teens[] = {
        3, 6, 6, 8, 8, 7, 7, 9, 8, 8,
    };
    static const int tens[] = {
        0, 0, 6, 6, 5, 5, 5, 7, 6, 6,
    };

    if (n == 1000) {
        return 11;  // "one thousand"
    }

    int count = 0;
    if (n >= 100) {
        count += ones[n / 100] + 7;  // "x hundred"
        n %= 100;
        if (n != 0) {
            count += 3;  // "and"
        }
    }

    if (n >= 20) {
        count += tens[n / 10] + ones[n % 10];
    } else if (n >= 10) {
        count += teens[n - 10];
    } else {
        count += ones[n];
    }

    return count;
}

int solve(const int limit) {
    int total = 0;
    for (int n = 1; n <= limit; ++n) {
        total += letters_for_number(n);
    }
    return total;
}

bool run_checkpoints() {
    if (solve(5) != 19) {
        std::cerr << "Checkpoint failed for limit=5" << '\n';
        return false;
    }
    if (letters_for_number(342) != 23) {
        std::cerr << "Checkpoint failed for n=342" << '\n';
        return false;
    }
    if (letters_for_number(115) != 20) {
        std::cerr << "Checkpoint failed for n=115" << '\n';
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
