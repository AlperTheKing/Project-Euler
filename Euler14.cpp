#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

struct Options {
    u64 limit = 1000000ULL;
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

    u64 parsed = 0ULL;
    for (const char c : tail) {
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

bool parse_arguments(const int argc, char** argv, Options& options) {
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

u32 collatz_length(const u64 start, std::unordered_map<u64, u32>& memo) {
    auto it = memo.find(start);
    if (it != memo.end()) {
        return it->second;
    }

    std::vector<u64> stack;
    u64 n = start;

    while (memo.find(n) == memo.end()) {
        stack.push_back(n);
        if ((n & 1ULL) == 0ULL) {
            n >>= 1ULL;
        } else {
            if (n > (std::numeric_limits<u64>::max() - 1ULL) / 3ULL) {
                throw std::overflow_error("Collatz overflow");
            }
            n = 3ULL * n + 1ULL;
        }
    }

    u32 length = memo[n];
    while (!stack.empty()) {
        const u64 x = stack.back();
        stack.pop_back();
        ++length;
        memo[x] = length;
    }

    return memo[start];
}

u64 solve(const u64 limit) {
    std::unordered_map<u64, u32> memo;
    memo.reserve(static_cast<std::size_t>(limit) * 2U);
    memo[1ULL] = 1U;

    u64 best_n = 1ULL;
    u32 best_len = 1U;

    for (u64 n = 2ULL; n < limit; ++n) {
        const u32 length = collatz_length(n, memo);
        if (length > best_len) {
            best_len = length;
            best_n = n;
        }
    }

    return best_n;
}

bool run_checkpoints() {
    std::unordered_map<u64, u32> memo;
    memo[1ULL] = 1U;
    if (collatz_length(13ULL, memo) != 10U) {
        std::cerr << "Checkpoint failed for start=13" << '\n';
        return false;
    }
    if (solve(20ULL) != 18ULL) {
        std::cerr << "Checkpoint failed for limit=20" << '\n';
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
