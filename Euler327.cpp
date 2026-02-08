#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u128 = unsigned __int128;

struct Options {
    int c_min = 3;
    int c_max = 40;
    int rooms = 30;
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
        if (parse_int_after_prefix(arg, "--c-min=", options.c_min) ||
            parse_int_after_prefix(arg, "--c-max=", options.c_max) ||
            parse_int_after_prefix(arg, "--rooms=", options.rooms)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.c_min >= 3 && options.c_max >= options.c_min && options.rooms >= 1;
}

u128 M(const int c, const int rooms) {
    u128 f = 1;  // M-1 helper.
    for (int r = 2; r <= rooms; ++r) {
        const u128 trips = (f + static_cast<u128>(c - 3)) / static_cast<u128>(c - 2);
        f = f + 2 * trips - 1;
    }
    return f + 1;
}

u128 solve(const int c_min, const int c_max, const int rooms) {
    u128 sum = 0;
    for (int c = c_min; c <= c_max; ++c) {
        sum += M(c, rooms);
    }
    return sum;
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
    if (M(3, 6) != 123 || M(4, 6) != 23) {
        std::cerr << "Checkpoint failed for M(3,6) and M(4,6) samples" << '\n';
        return false;
    }
    if (solve(3, 4, 6) != 146) {
        std::cerr << "Checkpoint failed for sum M(C,6), 3<=C<=4" << '\n';
        return false;
    }
    if (solve(3, 10, 10) != 10382) {
        std::cerr << "Checkpoint failed for sum M(C,10), 3<=C<=10" << '\n';
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
    std::cout << to_string_u128(solve(options.c_min, options.c_max, options.rooms)) << '\n';
    return 0;
}
