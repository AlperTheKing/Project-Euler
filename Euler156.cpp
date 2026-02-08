#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    i64 limit = 1000000000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 0;
}

i64 count_digit_upto(i64 n, const int d) {
    ++n;

    i64 result = 0;
    i64 p10 = 1;
    while (p10 < n) {
        const int cur = static_cast<int>((n / p10) % 10);
        result += (n / p10 / 10) * p10;

        if (cur == d) {
            result += n % (10 * p10) - static_cast<i64>(d) * p10;
        } else if (cur > d) {
            result += p10;
        }

        p10 *= 10;
    }

    return result;
}

void search_fixed_points(const int d, const i64 start, const i64 span, const i64 limit, i64& acc) {
    if (start > limit) {
        return;
    }

    const i64 right = std::min(limit, start + span);
    const i64 lo = count_digit_upto(start, d);
    const i64 hi = count_digit_upto(right, d);

    if (hi < start || lo > right) {
        return;
    }

    if (span == 1) {
        if (lo == start) {
            acc += start;
        }
        return;
    }

    const i64 step = span / 10;
    for (int k = 0; k < 10; ++k) {
        search_fixed_points(d, start + static_cast<i64>(k) * step, step, limit, acc);
    }
}

i64 solve(const i64 limit) {
    i64 root_span = 1;
    while (root_span < limit) {
        root_span *= 10;
    }

    i64 sum = 0;
    for (int d = 1; d <= 9; ++d) {
        search_fixed_points(d, 0, root_span, limit, sum);
    }
    return sum;
}

i64 brute(const i64 limit) {
    auto count_in_number = [](i64 x, const int d) {
        if (x == 0) {
            return 0;
        }
        int cnt = 0;
        while (x > 0) {
            cnt += static_cast<int>(x % 10 == d);
            x /= 10;
        }
        return cnt;
    };

    i64 ans = 0;
    for (int d = 1; d <= 9; ++d) {
        i64 running = 0;
        for (i64 x = 0; x <= limit; ++x) {
            running += count_in_number(x, d);
            if (running == x) {
                ans += x;
            }
        }
    }
    return ans;
}

bool run_checkpoints() {
    if (count_digit_upto(12, 1) != 5) {
        std::cerr << "Checkpoint failed for f(12,1)" << '\n';
        return false;
    }
    if (solve(100000) != brute(100000)) {
        std::cerr << "Checkpoint failed for brute cross-check at 100000" << '\n';
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
