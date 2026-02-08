#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 index = 1000000ULL;
    std::string digits = "0123456789";
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

bool parse_string_after_prefix(const std::string& arg,
                               const std::string& prefix,
                               std::string& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    value = tail;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--index=", options.index)) {
            continue;
        }
        if (parse_string_after_prefix(arg, "--digits=", options.digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.index >= 1ULL && !options.digits.empty();
}

u64 factorial(const std::size_t n) {
    u64 value = 1ULL;
    for (std::size_t k = 2; k <= n; ++k) {
        if (value > std::numeric_limits<u64>::max() / static_cast<u64>(k)) {
            throw std::overflow_error("factorial overflow");
        }
        value *= static_cast<u64>(k);
    }
    return value;
}

std::string nth_lexicographic_permutation(std::string digits, u64 index_one_based) {
    std::sort(digits.begin(), digits.end());
    const std::size_t n = digits.size();

    const u64 total = factorial(n);
    if (index_one_based < 1ULL || index_one_based > total) {
        throw std::runtime_error("Permutation index out of range");
    }

    u64 rank = index_one_based - 1ULL;
    std::string answer;
    answer.reserve(n);

    std::vector<char> pool(digits.begin(), digits.end());

    for (std::size_t remaining = n; remaining > 0; --remaining) {
        const u64 block = factorial(remaining - 1);
        const u64 pick = rank / block;
        rank %= block;

        answer.push_back(pool[static_cast<std::size_t>(pick)]);
        pool.erase(pool.begin() + static_cast<std::ptrdiff_t>(pick));
    }

    return answer;
}

std::string solve(const Options& options) {
    return nth_lexicographic_permutation(options.digits, options.index);
}

bool run_checkpoints() {
    if (nth_lexicographic_permutation("012", 1ULL) != "012") {
        std::cerr << "Checkpoint failed for first permutation" << '\n';
        return false;
    }
    if (nth_lexicographic_permutation("012", 6ULL) != "210") {
        std::cerr << "Checkpoint failed for last permutation" << '\n';
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

    try {
        std::cout << solve(options) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
