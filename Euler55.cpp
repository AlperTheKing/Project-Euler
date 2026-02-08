#include <algorithm>
#include <iostream>
#include <string>

namespace {

struct Options {
    int limit = 10000;
    int max_steps = 50;
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
        if (parse_int_after_prefix(arg, "--max-steps=", options.max_steps)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1 && options.max_steps >= 1;
}

bool is_palindrome(const std::string& s) {
    for (std::size_t i = 0; i < s.size() / 2U; ++i) {
        if (s[i] != s[s.size() - 1U - i]) {
            return false;
        }
    }
    return true;
}

std::string add_decimal_strings(const std::string& a, const std::string& b) {
    const std::size_t n = std::max(a.size(), b.size());
    std::string out;
    out.reserve(n + 1U);

    int carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        int sum = carry;
        if (i < a.size()) {
            sum += a[a.size() - 1U - i] - '0';
        }
        if (i < b.size()) {
            sum += b[b.size() - 1U - i] - '0';
        }
        out.push_back(static_cast<char>('0' + (sum % 10)));
        carry = sum / 10;
    }

    while (carry > 0) {
        out.push_back(static_cast<char>('0' + (carry % 10)));
        carry /= 10;
    }

    std::reverse(out.begin(), out.end());
    return out;
}

bool is_lychrel_candidate(const int n, const int max_steps) {
    std::string value = std::to_string(n);

    for (int step = 0; step < max_steps; ++step) {
        std::string rev = value;
        std::reverse(rev.begin(), rev.end());
        value = add_decimal_strings(value, rev);
        if (is_palindrome(value)) {
            return false;
        }
    }

    return true;
}

int solve(const int limit, const int max_steps) {
    int count = 0;
    for (int n = 1; n < limit; ++n) {
        if (is_lychrel_candidate(n, max_steps)) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (is_lychrel_candidate(47, 50)) {
        std::cerr << "Checkpoint failed for 47" << '\n';
        return false;
    }
    if (is_lychrel_candidate(349, 3)) {
        std::cerr << "Checkpoint failed for 349" << '\n';
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

    std::cout << solve(options.limit, options.max_steps) << '\n';
    return 0;
}
