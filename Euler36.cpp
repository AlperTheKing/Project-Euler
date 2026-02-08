#include <algorithm>
#include <iostream>
#include <string>

namespace {

struct Options {
    int limit = 1000000;
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
    for (const char c : tail) {
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

bool is_palindrome(const std::string& s) {
    for (std::size_t i = 0; i < s.size() / 2U; ++i) {
        if (s[i] != s[s.size() - 1U - i]) {
            return false;
        }
    }
    return true;
}

std::string to_binary(int value) {
    std::string out;
    while (value > 0) {
        out.push_back((value & 1) ? '1' : '0');
        value >>= 1;
    }
    if (out.empty()) {
        out = "0";
    }
    std::reverse(out.begin(), out.end());
    return out;
}

long long solve(const int limit) {
    long long total = 0;
    for (int n = 1; n < limit; ++n) {
        if (is_palindrome(std::to_string(n)) && is_palindrome(to_binary(n))) {
            total += n;
        }
    }
    return total;
}

bool run_checkpoints() {
    if (solve(1000) != 1772LL) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
        return false;
    }
    if (!is_palindrome("12321")) {
        std::cerr << "Checkpoint failed for palindrome helper" << '\n';
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
