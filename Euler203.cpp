#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int rows = 51;
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
        if (parse_int_after_prefix(arg, "--rows=", options.rows)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.rows >= 1;
}

bool is_squarefree(u64 v) {
    if (v == 0) {
        return false;
    }

    if ((v % 4ULL) == 0ULL) {
        return false;
    }
    for (u64 p = 3; p * p <= v; p += 2) {
        const u64 pp = p * p;
        if ((v % pp) == 0ULL) {
            return false;
        }
    }
    return true;
}

u64 solve(const int rows) {
    std::set<u64> values;

    for (int n = 0; n < rows; ++n) {
        u64 c = 1;
        for (int k = 0; k <= n; ++k) {
            if (k > 0) {
                c = c * static_cast<u64>(n - k + 1) / static_cast<u64>(k);
            }
            values.insert(c);
        }
    }

    u64 sum = 0;
    for (u64 v : values) {
        if (is_squarefree(v)) {
            sum += v;
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(8) != 105ULL) {
        std::cerr << "Checkpoint failed for first 8 rows" << '\n';
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

    std::cout << solve(options.rows) << '\n';
    return 0;
}
