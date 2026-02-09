#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    i64 perimeter_limit = 75'000'000;
    bool run_checkpoints = true;
};

struct Triple {
    i64 a;
    i64 b;
    i64 c;
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

Triple transform(const std::array<std::array<i64, 3>, 3>& m, const Triple& t) {
    const i64 x = m[0][0] * t.a + m[0][1] * t.b + m[0][2] * t.c;
    const i64 y = m[1][0] * t.a + m[1][1] * t.b + m[1][2] * t.c;
    const i64 z = m[2][0] * t.a + m[2][1] * t.b + m[2][2] * t.c;
    return {x, y, z};
}

u64 count_barely_obtuse(const i64 perimeter_limit) {
    static constexpr std::array<std::array<i64, 3>, 3> kM1{{
        {{1, -2, 2}},
        {{2, -1, 2}},
        {{2, -2, 3}},
    }};
    static constexpr std::array<std::array<i64, 3>, 3> kM2{{
        {{1, 2, 2}},
        {{2, 1, 2}},
        {{2, 2, 3}},
    }};
    static constexpr std::array<std::array<i64, 3>, 3> kM3{{
        {{-1, 2, 2}},
        {{-2, 1, 2}},
        {{-2, 2, 3}},
    }};

    std::vector<Triple> stack;
    stack.reserve(1 << 20);
    stack.push_back({2, 2, 3});

    u64 count = 0;

    while (!stack.empty()) {
        const Triple cur = stack.back();
        stack.pop_back();

        const i64 perimeter = cur.a + cur.b + cur.c;
        if (perimeter > perimeter_limit) {
            continue;
        }

        if (cur.a <= cur.b && cur.a + cur.b > cur.c) {
            ++count;
        }

        const Triple children[3] = {
            transform(kM1, cur),
            transform(kM2, cur),
            transform(kM3, cur),
        };

        for (const Triple& child : children) {
            if (child.a <= 0 || child.b <= 0 || child.c <= 0) {
                continue;
            }
            if (child.a + child.b + child.c <= perimeter_limit) {
                stack.push_back(child);
            }
        }
    }

    return count;
}

u64 brute_count(const int perimeter_limit) {
    u64 count = 0;
    for (int a = 1; a <= perimeter_limit / 3; ++a) {
        for (int b = a; b <= (perimeter_limit - a) / 2; ++b) {
            const i64 c2 = 1LL * a * a + 1LL * b * b + 1;
            i64 c = static_cast<i64>(std::sqrt(static_cast<long double>(c2)));
            while ((c + 1) * (c + 1) <= c2) {
                ++c;
            }
            while (c * c > c2) {
                --c;
            }
            if (c * c != c2) {
                continue;
            }
            if (a + b <= c) {
                continue;
            }
            if (a + b + c > perimeter_limit) {
                continue;
            }
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    for (int perimeter : {80, 120, 200, 500, 1000, 2000}) {
        const u64 slow = brute_count(perimeter);
        const u64 fast = count_barely_obtuse(perimeter);
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
    std::cout << count_barely_obtuse(options.perimeter_limit) << '\n';
    return 0;
}
