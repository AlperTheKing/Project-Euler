#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u128 = unsigned __int128;

struct Options {
    int index = 15;
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
        if (parse_int_after_prefix(arg, "--index=", options.index)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.index >= 1;
}

u128 nth_nugget(const int index) {
    const int need = 2 * index + 1;

    u128 f0 = 0;
    u128 f1 = 1;
    for (int i = 2; i <= need; ++i) {
        const u128 f2 = f0 + f1;
        f0 = f1;
        f1 = f2;
    }

    const u128 f_2n = f0;
    const u128 f_2n1 = f1;
    return f_2n * f_2n1;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool run_checkpoints() {
    if (to_string_u128(nth_nugget(3)) != "104") {
        std::cerr << "Checkpoint failed for index=3" << '\n';
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

    std::cout << to_string_u128(nth_nugget(options.index)) << '\n';
    return 0;
}
