#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int digits = 3;
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
        if (parse_int_after_prefix(arg, "--digits=", options.digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.digits > 0;
}

bool is_palindrome(i64 value) {
    std::string s = std::to_string(value);
    for (std::size_t i = 0; i < s.size() / 2U; ++i) {
        if (s[i] != s[s.size() - 1U - i]) {
            return false;
        }
    }
    return true;
}

i64 solve(const int digits) {
    i64 low = 1;
    for (int i = 1; i < digits; ++i) {
        low *= 10;
    }
    const i64 high = low * 10 - 1;

    i64 best = 0;
    for (i64 a = high; a >= low; --a) {
        if (a * high < best) {
            break;
        }

        for (i64 b = a; b >= low; --b) {
            const i64 prod = a * b;
            if (prod <= best) {
                break;
            }
            if (is_palindrome(prod)) {
                best = prod;
            }
        }
    }

    return best;
}

bool run_checkpoints() {
    if (solve(2) != 9009LL) {
        std::cerr << "Checkpoint failed for digits=2" << '\n';
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

    std::cout << solve(options.digits) << '\n';
    return 0;
}
