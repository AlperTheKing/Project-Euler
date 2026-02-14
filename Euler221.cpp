#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int target = 150000;
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
        if (parse_int_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 1;
}

u64 nth_alexandrian(const int target) {
    std::vector<u128> values;

    int p_limit = 2000;
    while (true) {
        values.clear();
        values.reserve(static_cast<std::size_t>(p_limit * 10));

        for (int p = 1; p <= p_limit; ++p) {
            const u64 n = static_cast<u64>(p) * static_cast<u64>(p) + 1ULL;
            const int r = static_cast<int>(std::sqrt(static_cast<long double>(n)));

            for (int d = 1; d <= r; ++d) {
                if (n % static_cast<u64>(d) != 0ULL) {
                    continue;
                }
                const u64 e = n / static_cast<u64>(d);
                const u128 pu = static_cast<u128>(static_cast<u64>(p));
                const u128 v = pu * (pu + static_cast<u128>(static_cast<u64>(d))) * (pu + static_cast<u128>(e));
                values.push_back(v);
            }
        }

        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());

        if (static_cast<int>(values.size()) >= target) {
            const u128 candidate = values[static_cast<std::size_t>(target - 1)];
            const u128 lb = static_cast<u128>(p_limit + 1) * static_cast<u128>(p_limit + 2) * static_cast<u128>(p_limit + 2);
            if (lb > candidate) {
                return static_cast<u64>(candidate);
            }
        }

        p_limit *= 2;
    }
}

bool run_checkpoints() {
    if (nth_alexandrian(6) != 630ULL) {
        std::cerr << "Checkpoint failed for 6th Alexandrian integer" << '\n';
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

    std::cout << nth_alexandrian(options.target) << '\n';
    return 0;
}
