#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int length = 20;
    int modulo = 1000000000;
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
        if (parse_int_after_prefix(arg, "--length=", options.length) ||
            parse_int_after_prefix(arg, "--modulo=", options.modulo)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.length >= 1 && options.modulo >= 1;
}

u64 solve(const int length, const int modulo) {
    const int max_sq_sum = length * 81;
    std::vector<std::vector<u64>> count(static_cast<std::size_t>(length + 1),
                                        std::vector<u64>(static_cast<std::size_t>(max_sq_sum + 1), 0));
    std::vector<std::vector<u64>> sum(static_cast<std::size_t>(length + 1),
                                      std::vector<u64>(static_cast<std::size_t>(max_sq_sum + 1), 0));
    std::vector<u64> pow10(static_cast<std::size_t>(length + 1), 1);
    for (int i = 1; i <= length; ++i) {
        pow10[static_cast<std::size_t>(i)] = (pow10[static_cast<std::size_t>(i - 1)] * 10ULL) %
                                             static_cast<u64>(modulo);
    }

    count[0][0] = 1;
    for (int len = 1; len <= length; ++len) {
        for (int digit = 0; digit <= 9; ++digit) {
            const int sq = digit * digit;
            for (int s = sq; s <= max_sq_sum; ++s) {
                const u64 cnt_prev = count[static_cast<std::size_t>(len - 1)][static_cast<std::size_t>(s - sq)];
                const u64 sum_prev = sum[static_cast<std::size_t>(len - 1)][static_cast<std::size_t>(s - sq)];
                if (cnt_prev == 0 && sum_prev == 0) {
                    continue;
                }
                count[static_cast<std::size_t>(len)][static_cast<std::size_t>(s)] =
                    (count[static_cast<std::size_t>(len)][static_cast<std::size_t>(s)] + cnt_prev) %
                    static_cast<u64>(modulo);

                const u64 add = (sum_prev +
                                 (pow10[static_cast<std::size_t>(len - 1)] *
                                  static_cast<u64>(digit) % static_cast<u64>(modulo) * cnt_prev) %
                                     static_cast<u64>(modulo)) %
                                static_cast<u64>(modulo);
                sum[static_cast<std::size_t>(len)][static_cast<std::size_t>(s)] =
                    (sum[static_cast<std::size_t>(len)][static_cast<std::size_t>(s)] + add) %
                    static_cast<u64>(modulo);
            }
        }
    }

    u64 answer = 0;
    for (int r = 1; r * r <= max_sq_sum; ++r) {
        answer = (answer + sum[static_cast<std::size_t>(length)][static_cast<std::size_t>(r * r)]) %
                 static_cast<u64>(modulo);
    }
    return answer;
}

u64 brute_small(const int length, const int modulo) {
    u64 answer = 0;
    int total_values = 1;
    for (int i = 0; i < length; ++i) {
        total_values *= 10;
    }

    for (int x = 1; x < total_values; ++x) {
        int y = x;
        int sq_sum = 0;
        for (int i = 0; i < length; ++i) {
            const int d = y % 10;
            sq_sum += d * d;
            y /= 10;
        }
        const int r = static_cast<int>(std::sqrt(static_cast<double>(sq_sum)));
        if (r * r == sq_sum) {
            answer = (answer + static_cast<u64>(x)) % static_cast<u64>(modulo);
        }
    }
    return answer;
}

bool run_checkpoints() {
    if (solve(5, 1000000000) != brute_small(5, 1000000000)) {
        std::cerr << "Checkpoint failed for length 5 brute cross-check" << '\n';
        return false;
    }
    if (solve(6, 1000000000) != brute_small(6, 1000000000)) {
        std::cerr << "Checkpoint failed for length 6 brute cross-check" << '\n';
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
    std::cout << solve(options.length, options.modulo) << '\n';
    return 0;
}
