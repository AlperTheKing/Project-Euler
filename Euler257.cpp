#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    i64 perimeter_limit = 100'000'000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_int_after_prefix(arg, "--perimeter=", options.perimeter_limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.perimeter_limit >= 3;
}

i64 floor_sqrt(i64 x) {
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

i64 ceil_sqrt(i64 x) {
    const i64 r = floor_sqrt(x);
    return (r * r == x) ? r : (r + 1);
}

u64 count_valid_triangles(const i64 perimeter_limit) {
    u64 count = static_cast<u64>(perimeter_limit / 3);  // k = 4 (equilateral)

    const i64 max_var = floor_sqrt(perimeter_limit) + 2;

    // k = 2, decomposition (b-a)(c-a) = 2a^2
    // canonical case A: d1 = g*x^2, d2 = 2g*y^2, gcd(x,y)=1, x odd
    for (i64 y = 1; y <= max_var; ++y) {
        const i64 x_max = floor_sqrt(2 * y * y);
        for (i64 x = y + 1; x <= x_max; ++x) {
            if (std::gcd(x, y) != 1 || (x & 1LL) == 0) {
                continue;
            }
            const i64 base = (x + y) * (x + 2 * y);
            if (base > perimeter_limit) {
                continue;
            }
            count += static_cast<u64>(perimeter_limit / base);
        }
    }

    // canonical case B: d1 = 2g*x^2, d2 = g*y^2, gcd(x,y)=1, y odd
    for (i64 x = 1; x <= max_var; ++x) {
        const i64 y_min = ceil_sqrt(2 * x * x);
        const i64 y_max = 2 * x - 1;
        for (i64 y = y_min; y <= y_max; ++y) {
            if (std::gcd(x, y) != 1 || (y & 1LL) == 0) {
                continue;
            }
            const i64 base = (x + y) * (2 * x + y);
            if (base > perimeter_limit) {
                continue;
            }
            count += static_cast<u64>(perimeter_limit / base);
        }
    }

    // k = 3, decomposition (2b-a)(2c-a) = 3a^2
    // canonical case A: u = g*x^2, v = 3g*y^2, gcd(x,y)=1, x % 3 != 0
    for (i64 y = 1; y <= max_var; ++y) {
        const i64 x_max = floor_sqrt(3 * y * y);
        for (i64 x = y + 1; x <= x_max; ++x) {
            if (std::gcd(x, y) != 1 || x % 3 == 0) {
                continue;
            }

            const i64 base_num = (x + y) * (x + 3 * y);  // perimeter = g * base_num / 2
            if (base_num > 2 * perimeter_limit) {
                continue;
            }

            const i64 g_max = (2 * perimeter_limit) / base_num;
            const bool odd_g_valid = ((x * (x + y)) % 2 == 0) && ((y * (x + 3 * y)) % 2 == 0);
            count += odd_g_valid ? static_cast<u64>(g_max) : static_cast<u64>(g_max / 2);
        }
    }

    // canonical case B: u = 3g*x^2, v = g*y^2, gcd(x,y)=1, y % 3 != 0
    for (i64 x = 1; x <= max_var; ++x) {
        const i64 y_min = ceil_sqrt(3 * x * x);
        const i64 y_max = 3 * x - 1;
        for (i64 y = y_min; y <= y_max; ++y) {
            if (std::gcd(x, y) != 1 || y % 3 == 0) {
                continue;
            }

            const i64 base_num = (x + y) * (3 * x + y);  // perimeter = g * base_num / 2
            if (base_num > 2 * perimeter_limit) {
                continue;
            }

            const i64 g_max = (2 * perimeter_limit) / base_num;
            const bool odd_g_valid = ((x * (y + 3 * x)) % 2 == 0) && ((y * (x + y)) % 2 == 0);
            count += odd_g_valid ? static_cast<u64>(g_max) : static_cast<u64>(g_max / 2);
        }
    }

    return count;
}

u64 brute_count(const int perimeter_limit) {
    u64 count = 0;
    for (int a = 1; a <= perimeter_limit / 3; ++a) {
        for (int b = a; b <= (perimeter_limit - a) / 2; ++b) {
            for (int c = b; a + b + c <= perimeter_limit; ++c) {
                if (a + b <= c) {
                    continue;
                }
                const i64 num = static_cast<i64>(a + b) * static_cast<i64>(a + c);
                const i64 den = static_cast<i64>(b) * static_cast<i64>(c);
                if (num % den == 0) {
                    ++count;
                }
            }
        }
    }
    return count;
}

bool run_checkpoints() {
    for (int perimeter : {40, 60, 100, 150, 200, 300}) {
        const u64 slow = brute_count(perimeter);
        const u64 fast = count_valid_triangles(perimeter);
        if (slow != fast) {
            std::cerr << "Checkpoint failed at perimeter=" << perimeter
                      << ": brute=" << slow << ", fast=" << fast << '\n';
            return false;
        }
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
    std::cout << count_valid_triangles(options.perimeter_limit) << '\n';
    return 0;
}
