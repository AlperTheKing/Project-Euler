#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

struct Options {
    u64 n = 1'000'000'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n == 0ULL) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

std::vector<int> build_rows(const u64 n) {
    std::vector<int> rows;
    u64 p3 = 1ULL;
    while (p3 <= n) {
        const u64 bound = n / p3;
        const int max_a = 63 - __builtin_clzll(bound);
        rows.push_back(max_a + 1);
        if (p3 > n / 3ULL) {
            break;
        }
        p3 *= 3ULL;
    }
    return rows;
}

cpp_int count_exact(const u64 n) {
    const std::vector<int> rows = build_rows(n);
    const int max_width = rows.empty() ? 0 : rows.front();
    int cell_count = 0;
    for (const int w : rows) {
        cell_count += w;
    }

    std::vector<int> col_heights(static_cast<std::size_t>(max_width), 0);
    for (int c = 0; c < max_width; ++c) {
        int h = 0;
        for (const int w : rows) {
            if (w > c) {
                ++h;
            }
        }
        col_heights[static_cast<std::size_t>(c)] = h;
    }

    cpp_int numerator = 1;
    for (int x = 2; x <= cell_count; ++x) {
        numerator *= x;
    }

    cpp_int denominator = 1;
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        const int w = rows[static_cast<std::size_t>(r)];
        for (int c = 0; c < w; ++c) {
            const int right = w - c - 1;
            const int below = col_heights[static_cast<std::size_t>(c)] - r - 1;
            const int hook = right + below + 1;
            denominator *= hook;
        }
    }

    if ((numerator % denominator) != 0) {
        std::cerr << "Internal error: non-integral hook-length ratio.\n";
        return 0;
    }
    return numerator / denominator;
}

std::string format_scientific_10(const cpp_int& value) {
    std::string digits = value.convert_to<std::string>();
    long long exponent = static_cast<long long>(digits.size()) - 1LL;

    if (digits.size() < 12U) {
        digits.append(12U - digits.size(), '0');
    }

    u64 leading = 0ULL;
    for (int i = 0; i < 11; ++i) {
        leading = leading * 10ULL + static_cast<u64>(digits[static_cast<std::size_t>(i)] - '0');
    }
    const int round_digit = digits[11] - '0';
    if (round_digit >= 5) {
        ++leading;
    }
    if (leading == 100'000'000'000ULL) {
        leading = 10'000'000'000ULL;
        ++exponent;
    }

    const u64 integer_part = leading / 10'000'000'000ULL;
    const u64 fractional_part = leading % 10'000'000'000ULL;

    std::ostringstream out;
    out << integer_part << '.'
        << std::setw(10) << std::setfill('0') << fractional_part
        << 'e' << exponent;
    return out.str();
}

std::string solve_formatted(const u64 n) {
    return format_scientific_10(count_exact(n));
}

bool run_checkpoints() {
    if (count_exact(6ULL) != cpp_int(5)) {
        std::cerr << "Checkpoint failed: F(6)\n";
        return false;
    }
    if (count_exact(8ULL) != cpp_int(9)) {
        std::cerr << "Checkpoint failed: F(8)\n";
        return false;
    }
    if (count_exact(20ULL) != cpp_int(450)) {
        std::cerr << "Checkpoint failed: F(20)\n";
        return false;
    }
    if (solve_formatted(1000ULL) != "8.8521816557e21") {
        std::cerr << "Checkpoint failed: F(1000)\n";
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
        return 1;
    }

    std::cout << solve_formatted(options.n) << '\n';
    return 0;
}
