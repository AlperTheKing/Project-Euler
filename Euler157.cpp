#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = long long;

struct Options {
    int max_n = 9;
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
        if (parse_int_after_prefix(arg, "--max-n=", options.max_n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.max_n >= 1;
}

i64 ipow(const int base, int exp) {
    i64 out = 1;
    while (exp > 0) {
        out *= base;
        --exp;
    }
    return out;
}

int divisor_count(i64 x) {
    int count = 1;
    for (i64 p = 2; p * p <= x; ++p) {
        int e = 1;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        count *= e;
    }
    if (x != 1) {
        count *= 2;
    }
    return count;
}

i64 count_for_n(const int n) {
    i64 total = 0;
    for (int a2 = 0; a2 <= 0; ++a2) {
        for (int b2 = 0; a2 + b2 <= n; ++b2) {
            if (a2 != 0 && b2 != 0) {
                continue;
            }
            for (int a5 = 0; a5 <= (b2 ? n : 0); ++a5) {
                for (int b5 = 0; a5 + b5 <= n; ++b5) {
                    if (a5 != 0 && b5 != 0) {
                        continue;
                    }

                    const i64 left = ipow(2, a2) * ipow(5, a5);
                    const i64 right = ipow(2, b2) * ipow(5, b5);
                    const i64 scale = ipow(2, n - a2 - b2) * ipow(5, n - a5 - b5);
                    total += divisor_count((left + right) * scale);
                }
            }
        }
    }
    return total;
}

i64 solve(const int max_n) {
    i64 total = 0;
    for (int n = 1; n <= max_n; ++n) {
        total += count_for_n(n);
    }
    return total;
}

bool run_checkpoints() {
    if (count_for_n(1) != 20) {
        std::cerr << "Checkpoint failed for n=1" << '\n';
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

    std::cout << solve(options.max_n) << '\n';
    return 0;
}
