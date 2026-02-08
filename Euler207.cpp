#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 denominator = 12345ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--denominator=", options.denominator)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.denominator >= 1ULL;
}

u64 solve(const u64 denominator) {
    u64 p = 1;
    u64 next_pow2 = 4;

    for (u64 m = 2;; ++m) {
        if (m == next_pow2) {
            ++p;
            next_pow2 <<= 1;
        }

        if (static_cast<u128>(p) * static_cast<u128>(denominator) < static_cast<u128>(m - 1)) {
            return m * (m - 1);
        }
    }
}

bool run_checkpoints() {
    if (solve(3ULL) != 110ULL) {
        std::cerr << "Checkpoint failed for threshold 1/3" << '\n';
        return false;
    }
    if (solve(10ULL) != 2652ULL) {
        std::cerr << "Checkpoint failed for threshold 1/10" << '\n';
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

    std::cout << solve(options.denominator) << '\n';
    return 0;
}
