#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 limit = 1'000'000'000'000ULL;
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
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2ULL;
}

std::vector<u64> strong_repunits_below(const u64 limit) {
    std::vector<u64> values;
    if (limit > 1ULL) {
        values.push_back(1ULL);
    }

    for (u64 base = 2ULL;; ++base) {
        if (base > (std::numeric_limits<u64>::max() - 1ULL) / base) {
            break;
        }
        const u64 square = base * base;
        if (1ULL + base + square >= limit) {
            break;
        }

        u64 repunit = 1ULL + base + square;  // Length 3.
        while (repunit < limit) {
            values.push_back(repunit);
            if (repunit > (limit - 1ULL) / base) {
                break;
            }
            repunit = repunit * base + 1ULL;
        }
    }

    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

u64 solve(const u64 limit) {
    const std::vector<u64> values = strong_repunits_below(limit);
    u64 sum = 0ULL;
    for (const u64 v : values) {
        sum += v;
    }
    return sum;
}

bool run_checkpoints() {
    const std::vector<u64> sample = strong_repunits_below(50ULL);
    const std::vector<u64> expected = {1ULL, 7ULL, 13ULL, 15ULL, 21ULL, 31ULL, 40ULL, 43ULL};
    if (sample != expected) {
        std::cerr << "Checkpoint failed for limit=50 strong repunits list" << '\n';
        return false;
    }
    if (solve(1000ULL) != 15864ULL) {
        std::cerr << "Checkpoint failed for limit=1000 sample sum" << '\n';
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
