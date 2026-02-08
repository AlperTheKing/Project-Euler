#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u128 = unsigned __int128;

struct Options {
    int max_digits = 16;
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
        if (parse_int_after_prefix(arg, "--max-digits=", options.max_digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.max_digits >= 1;
}

u128 pow_u128(const u128 base, int exp) {
    u128 result = 1;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

u128 count_for_length(const int len) {
    const u128 total = 15U * pow_u128(16U, len - 1);

    const u128 miss0 = pow_u128(15U, len);
    const u128 miss1 = 14U * pow_u128(15U, len - 1);
    const u128 missA = 14U * pow_u128(15U, len - 1);

    const u128 miss01 = pow_u128(14U, len);
    const u128 miss0A = pow_u128(14U, len);
    const u128 miss1A = 13U * pow_u128(14U, len - 1);

    const u128 miss01A = pow_u128(13U, len);

    return total - (miss0 + miss1 + missA) + (miss01 + miss0A + miss1A) - miss01A;
}

u128 solve_value(const int max_digits) {
    u128 total = 0;
    for (int len = 1; len <= max_digits; ++len) {
        total += count_for_length(len);
    }
    return total;
}

std::string to_hex(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value & 0xFULL);
        if (digit < 10) {
            out.push_back(static_cast<char>('0' + digit));
        } else {
            out.push_back(static_cast<char>('A' + (digit - 10)));
        }
        value >>= 4;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool run_checkpoints() {
    if (count_for_length(3) != 4U) {
        std::cerr << "Checkpoint failed for 3-digit hexadecimal numbers" << '\n';
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

    std::cout << to_hex(solve_value(options.max_digits)) << '\n';
    return 0;
}
