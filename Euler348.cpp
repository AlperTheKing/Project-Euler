#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int target_count = 5;
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
        if (parse_int_after_prefix(arg, "--target-count=", options.target_count)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target_count >= 1;
}

u64 isqrt_u64(const u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / (x + 1ULL)) {
        ++x;
    }
    while (x > 0ULL && x > n / x) {
        --x;
    }
    return x;
}

bool is_palindrome(u64 value) {
    u64 reversed = 0ULL;
    u64 x = value;
    while (x > 0ULL) {
        reversed = reversed * 10ULL + (x % 10ULL);
        x /= 10ULL;
    }
    return reversed == value;
}

u64 pow10_u64(int exp) {
    u64 value = 1ULL;
    for (int i = 0; i < exp; ++i) {
        value *= 10ULL;
    }
    return value;
}

u64 make_palindrome(const u64 prefix, const bool odd_length) {
    u64 tail = odd_length ? prefix / 10ULL : prefix;
    u64 palindrome = prefix;
    while (tail > 0ULL) {
        palindrome = palindrome * 10ULL + (tail % 10ULL);
        tail /= 10ULL;
    }
    return palindrome;
}

void ensure_cubes_up_to(const u64 n, std::vector<u64>& cubes, u64& next_base) {
    while (true) {
        const u64 cube = next_base * next_base * next_base;
        if (cube + 4ULL > n) {
            break;
        }
        cubes.push_back(cube);
        ++next_base;
    }
}

int count_representations(const u64 n, const std::vector<u64>& cubes) {
    int count = 0;
    for (const u64 cube : cubes) {
        if (cube + 4ULL > n) {
            break;
        }
        const u64 remaining = n - cube;
        const u64 root = isqrt_u64(remaining);
        if (root >= 2ULL && root * root == remaining) {
            ++count;
            if (count > 4) {
                break;
            }
        }
    }
    return count;
}

std::vector<u64> brute_palindromes_with_exact_count(const u64 limit, const int exact_count) {
    std::unordered_map<u64, int> counts;
    for (u64 a = 2ULL; a * a <= limit; ++a) {
        const u64 square = a * a;
        for (u64 b = 2ULL;; ++b) {
            const u64 cube = b * b * b;
            if (square + cube > limit) {
                break;
            }
            const u64 value = square + cube;
            if (is_palindrome(value)) {
                ++counts[value];
            }
        }
    }
    std::vector<u64> result;
    for (const auto& row : counts) {
        if (row.second == exact_count) {
            result.push_back(row.first);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<u64> fast_palindromes_with_exact_count(const u64 limit, const int exact_count) {
    std::vector<u64> cubes;
    u64 next_base = 2ULL;
    ensure_cubes_up_to(limit, cubes, next_base);

    std::vector<u64> result;
    for (u64 value = 12ULL; value <= limit; ++value) {
        if (!is_palindrome(value)) {
            continue;
        }
        if (count_representations(value, cubes) == exact_count) {
            result.push_back(value);
        }
    }
    return result;
}

u64 solve(const int target_count) {
    std::vector<u64> cubes;
    u64 next_base = 2ULL;
    std::vector<u64> found;

    for (int length = 1; static_cast<int>(found.size()) < target_count; ++length) {
        const bool odd_length = (length % 2) == 1;
        const int prefix_digits = (length + 1) / 2;
        const u64 begin = (length == 1) ? 1ULL : pow10_u64(prefix_digits - 1);
        const u64 end = pow10_u64(prefix_digits) - 1ULL;

        for (u64 prefix = begin; prefix <= end; ++prefix) {
            const u64 palindrome = make_palindrome(prefix, odd_length);
            if (palindrome < 12ULL) {
                continue;
            }
            ensure_cubes_up_to(palindrome, cubes, next_base);
            if (count_representations(palindrome, cubes) == 4) {
                found.push_back(palindrome);
                if (static_cast<int>(found.size()) == target_count) {
                    break;
                }
            }
        }
    }

    u64 sum = 0ULL;
    for (const u64 v : found) {
        sum += v;
    }
    return sum;
}

bool run_checkpoints() {
    std::vector<u64> cubes;
    u64 next_base = 2ULL;
    ensure_cubes_up_to(5'229'225ULL, cubes, next_base);
    if (count_representations(5'229'225ULL, cubes) != 4) {
        std::cerr << "Checkpoint failed for statement palindrome 5229225" << '\n';
        return false;
    }

    const u64 small_limit = 200'000ULL;
    const std::vector<u64> brute = brute_palindromes_with_exact_count(small_limit, 4);
    const std::vector<u64> fast = fast_palindromes_with_exact_count(small_limit, 4);
    if (brute != fast) {
        std::cerr << "Checkpoint failed for brute-force cross-check at limit=200000" << '\n';
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
    std::cout << solve(options.target_count) << '\n';
    return 0;
}
