#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int digits = 1000;
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
        if (parse_int_after_prefix(arg, "--digits=", options.digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.digits >= 1;
}

std::vector<int> add_bigints(const std::vector<int>& a, const std::vector<int>& b) {
    const std::size_t len = std::max(a.size(), b.size());
    std::vector<int> result;
    result.reserve(len + 1U);

    int carry = 0;
    for (std::size_t i = 0; i < len; ++i) {
        int value = carry;
        if (i < a.size()) {
            value += a[i];
        }
        if (i < b.size()) {
            value += b[i];
        }

        result.push_back(value % 10);
        carry = value / 10;
    }

    while (carry > 0) {
        result.push_back(carry % 10);
        carry /= 10;
    }

    return result;
}

int solve(const int digits) {
    if (digits == 1) {
        return 1;
    }

    std::vector<int> f_prev(1, 1);  // F1
    std::vector<int> f_curr(1, 1);  // F2

    int index = 2;
    while (static_cast<int>(f_curr.size()) < digits) {
        std::vector<int> f_next = add_bigints(f_prev, f_curr);
        f_prev.swap(f_curr);
        f_curr.swap(f_next);
        ++index;
    }

    return index;
}

bool run_checkpoints() {
    if (solve(2) != 7) {
        std::cerr << "Checkpoint failed for digits=2" << '\n';
        return false;
    }
    if (solve(3) != 12) {
        std::cerr << "Checkpoint failed for digits=3" << '\n';
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

    std::cout << solve(options.digits) << '\n';
    return 0;
}
