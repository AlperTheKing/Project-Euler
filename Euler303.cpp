#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int limit = 10000;
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

std::string smallest_multiple_with_digits_leq_2(const int n) {
    if (n == 1) {
        return "1";
    }

    std::vector<int> parent(static_cast<std::size_t>(n), -1);
    std::vector<int> parent_digit(static_cast<std::size_t>(n), -1);
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(n), 0U);
    std::queue<int> q;

    for (int first_digit = 1; first_digit <= 2; ++first_digit) {
        const int rem = first_digit % n;
        if (seen[static_cast<std::size_t>(rem)] != 0U) {
            continue;
        }
        seen[static_cast<std::size_t>(rem)] = 1U;
        parent[static_cast<std::size_t>(rem)] = -2;  // root marker
        parent_digit[static_cast<std::size_t>(rem)] = first_digit;
        q.push(rem);
    }

    int end_rem = -1;
    while (!q.empty()) {
        const int rem = q.front();
        q.pop();
        if (rem == 0) {
            end_rem = rem;
            break;
        }

        for (int digit = 0; digit <= 2; ++digit) {
            const int next = (rem * 10 + digit) % n;
            if (seen[static_cast<std::size_t>(next)] != 0U) {
                continue;
            }
            seen[static_cast<std::size_t>(next)] = 1U;
            parent[static_cast<std::size_t>(next)] = rem;
            parent_digit[static_cast<std::size_t>(next)] = digit;
            q.push(next);
        }
    }

    std::string digits;
    int cur = end_rem;
    while (cur >= 0) {
        digits.push_back(static_cast<char>('0' + parent_digit[static_cast<std::size_t>(cur)]));
        cur = parent[static_cast<std::size_t>(cur)];
    }
    std::reverse(digits.begin(), digits.end());
    return digits;
}

cpp_int solve(const int limit) {
    cpp_int total = 0;
    for (int n = 1; n <= limit; ++n) {
        const std::string multiple = smallest_multiple_with_digits_leq_2(n);
        cpp_int value = 0;
        for (char c : multiple) {
            value *= 10;
            value += static_cast<int>(c - '0');
        }
        total += value / n;
    }
    return total;
}

bool run_checkpoints() {
    if (smallest_multiple_with_digits_leq_2(2) != "2") {
        std::cerr << "Checkpoint failed for f(2)=2" << '\n';
        return false;
    }
    if (smallest_multiple_with_digits_leq_2(3) != "12") {
        std::cerr << "Checkpoint failed for f(3)=12" << '\n';
        return false;
    }
    if (smallest_multiple_with_digits_leq_2(7) != "21") {
        std::cerr << "Checkpoint failed for f(7)=21" << '\n';
        return false;
    }
    if (smallest_multiple_with_digits_leq_2(42) != "210") {
        std::cerr << "Checkpoint failed for f(42)=210" << '\n';
        return false;
    }
    if (solve(100) != cpp_int("11363107")) {
        std::cerr << "Checkpoint failed for limit=100 sample sum" << '\n';
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
