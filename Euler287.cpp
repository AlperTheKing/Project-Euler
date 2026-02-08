#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int n = 24;
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.n <= 30;
}

u64 encode_length_fast(const int n) {
    const int size = 1 << n;
    const int c = 1 << (n - 1);
    const u64 r2 = 1ULL << (2 * n - 2);

    const auto rec = [&](auto&& self, const int x0, const int y0, const int side) -> u64 {
        const int x1 = x0 + side - 1;
        const int y1 = y0 + side - 1;

        const int nearest_x = std::clamp(c, x0, x1);
        const int nearest_y = std::clamp(c, y0, y1);
        const long long dx_min = static_cast<long long>(nearest_x) - c;
        const long long dy_min = static_cast<long long>(nearest_y) - c;
        const u64 min_d2 = static_cast<u64>(dx_min * dx_min + dy_min * dy_min);

        const long long dx_max = std::max(std::llabs(static_cast<long long>(x0) - c),
                                          std::llabs(static_cast<long long>(x1) - c));
        const long long dy_max = std::max(std::llabs(static_cast<long long>(y0) - c),
                                          std::llabs(static_cast<long long>(y1) - c));
        const u64 max_d2 = static_cast<u64>(dx_max * dx_max + dy_max * dy_max);

        if (max_d2 <= r2 || min_d2 > r2) {
            return 2ULL;
        }

        const int half = side / 2;
        return 1ULL + self(self, x0, y0 + half, half) + self(self, x0 + half, y0 + half, half) +
               self(self, x0, y0, half) + self(self, x0 + half, y0, half);
    };

    return rec(rec, 0, 0, size);
}

u64 encode_length_bruteforce(const int n) {
    const int size = 1 << n;
    const int c = 1 << (n - 1);
    const u64 r2 = 1ULL << (2 * n - 2);

    std::vector<std::vector<int>> black(static_cast<std::size_t>(size),
                                        std::vector<int>(static_cast<std::size_t>(size), 0));
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            const long long dx = static_cast<long long>(x) - c;
            const long long dy = static_cast<long long>(y) - c;
            black[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)] =
                (static_cast<u64>(dx * dx + dy * dy) <= r2) ? 1 : 0;
        }
    }

    std::vector<std::vector<int>> ps(static_cast<std::size_t>(size + 1),
                                     std::vector<int>(static_cast<std::size_t>(size + 1), 0));
    for (int y = 0; y < size; ++y) {
        int row_sum = 0;
        for (int x = 0; x < size; ++x) {
            row_sum += black[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)];
            ps[static_cast<std::size_t>(y + 1)][static_cast<std::size_t>(x + 1)] =
                ps[static_cast<std::size_t>(y)][static_cast<std::size_t>(x + 1)] + row_sum;
        }
    }

    const auto sum_black = [&](const int x0, const int y0, const int side) {
        const int x1 = x0 + side;
        const int y1 = y0 + side;
        return ps[static_cast<std::size_t>(y1)][static_cast<std::size_t>(x1)] -
               ps[static_cast<std::size_t>(y0)][static_cast<std::size_t>(x1)] -
               ps[static_cast<std::size_t>(y1)][static_cast<std::size_t>(x0)] +
               ps[static_cast<std::size_t>(y0)][static_cast<std::size_t>(x0)];
    };

    const auto rec = [&](auto&& self, const int x0, const int y0, const int side) -> u64 {
        const int total = side * side;
        const int blacks = sum_black(x0, y0, side);
        if (blacks == 0 || blacks == total) {
            return 2ULL;
        }
        const int half = side / 2;
        return 1ULL + self(self, x0, y0 + half, half) + self(self, x0 + half, y0 + half, half) +
               self(self, x0, y0, half) + self(self, x0 + half, y0, half);
    };

    return rec(rec, 0, 0, size);
}

bool run_checkpoints() {
    for (int n = 1; n <= 5; ++n) {
        if (encode_length_fast(n) != encode_length_bruteforce(n)) {
            std::cerr << "Checkpoint failed for N=" << n << '\n';
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
    std::cout << encode_length_fast(options.n) << '\n';
    return 0;
}
