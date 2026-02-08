#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 100000000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

u64 sigma_prefix(const int n, std::unordered_map<int, u64>& memo) {
    const auto it = memo.find(n);
    if (it != memo.end()) {
        return it->second;
    }

    u64 sum = 0;
    int l = 1;
    while (l <= n) {
        const int q = n / l;
        const int r = n / q;
        const u64 interval_sum = (static_cast<u64>(l) + static_cast<u64>(r)) *
                                 static_cast<u64>(r - l + 1) / 2ULL;
        sum += static_cast<u64>(q) * interval_sum;
        l = r + 1;
    }
    memo.emplace(n, sum);
    return sum;
}

u64 solve(const int limit) {
    std::unordered_map<int, u64> sigma_cache;
    sigma_cache.reserve(32768);

    u64 answer = sigma_prefix(limit, sigma_cache);
    for (int a = 1; static_cast<u64>(a) * static_cast<u64>(a) <= static_cast<u64>(limit); ++a) {
        const int aa = a * a;
        const int bmax = static_cast<int>(std::sqrt(static_cast<double>(limit - aa)));
        for (int b = 1; b <= bmax; ++b) {
            if (std::gcd(a, b) != 1) {
                continue;
            }
            const int norm = aa + b * b;
            answer += static_cast<u64>(2 * a) * sigma_prefix(limit / norm, sigma_cache);
        }
    }
    return answer;
}

bool run_checkpoints() {
    if (solve(5) != 35ULL) {
        std::cerr << "Checkpoint failed for limit 5" << '\n';
        return false;
    }
    if (solve(100000) != 17924657155ULL) {
        std::cerr << "Checkpoint failed for limit 100000" << '\n';
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
