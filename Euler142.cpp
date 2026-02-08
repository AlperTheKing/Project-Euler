#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = long long;

struct Options {
    int search_limit = 5000;
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
        if (parse_int_after_prefix(arg, "--search-limit=", options.search_limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.search_limit >= 4;
}

bool is_square(const i64 x) {
    if (x < 0) {
        return false;
    }
    const i64 root = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    return root * root == x || (root + 1) * (root + 1) == x;
}

i64 solve(const int search_limit) {
    for (i64 i = 4; i <= search_limit; ++i) {
        const i64 a = i * i;
        for (i64 j = 3; j < i; ++j) {
            const i64 c = j * j;
            const i64 f = a - c;
            if (f <= 0 || !is_square(f)) {
                continue;
            }

            const i64 kstart = (j & 1LL) ? 1 : 2;
            for (i64 k = kstart; k < j; k += 2) {
                const i64 d = k * k;
                const i64 e = a - d;
                const i64 b = c - e;

                if (b <= 0 || e <= 0 || !is_square(b) || !is_square(e)) {
                    continue;
                }

                const i64 x = (d + c) / 2;
                const i64 y = (e + f) / 2;
                const i64 z = (c - d) / 2;
                if (!(x > y && y > z && z > 0)) {
                    continue;
                }

                if (is_square(x + y) && is_square(x - y) &&
                    is_square(x + z) && is_square(x - z) &&
                    is_square(y + z) && is_square(y - z)) {
                    return x + y + z;
                }
            }
        }
    }

    return 0;
}

bool run_checkpoints() {
    if (solve(5000) != 1006193) {
        std::cerr << "Checkpoint failed for search limit 5000" << '\n';
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

    std::cout << solve(options.search_limit) << '\n';
    return 0;
}
