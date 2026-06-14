#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 MOD = 1'003'443'221LL;

struct Interval {
    int left;
    int right;
};

std::vector<int> load_csv(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) {
        throw std::runtime_error("cannot open input file");
    }

    std::string data;
    std::getline(fin, data);
    std::stringstream ss(data);
    std::string token;
    std::vector<int> values;

    while (std::getline(ss, token, ',')) {
        if (!token.empty()) {
            values.push_back(std::stoi(token));
        }
    }

    return values;
}

std::vector<Interval> build_intervals(const std::vector<int>& values) {
    std::unordered_map<int, std::vector<int>> positions;
    positions.reserve(values.size() / 2U + 1U);

    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        positions[values[static_cast<std::size_t>(i)]].push_back(i);
    }

    std::vector<Interval> intervals;
    intervals.reserve(positions.size());
    for (const auto& entry : positions) {
        assert(entry.second.size() == 2U);
        intervals.push_back({entry.second[0], entry.second[1]});
    }

    std::sort(intervals.begin(), intervals.end(), [](const Interval& a, const Interval& b) {
        return a.left < b.left;
    });

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        assert(intervals[i - 1].left < intervals[i].left);
    }

    return intervals;
}

i64 add_mod(const i64 a, const i64 b, const i64 mod) {
    const i64 s = a + b;
    return s >= mod ? s - mod : s;
}

i64 mul_mod(const i64 a, const i64 b, const i64 mod) {
    return static_cast<i64>((static_cast<__int128>(a) * b) % mod);
}

i64 connectivity_number(const std::vector<int>& values, const i64 mod) {
    const std::vector<Interval> intervals = build_intervals(values);
    const int n = static_cast<int>(intervals.size());

    std::vector<int> left(static_cast<std::size_t>(n));
    std::vector<int> right(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        left[static_cast<std::size_t>(i)] = intervals[static_cast<std::size_t>(i)].left;
        right[static_cast<std::size_t>(i)] = intervals[static_cast<std::size_t>(i)].right;
    }

    std::vector<int> next(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        next[static_cast<std::size_t>(i)] = static_cast<int>(
            std::upper_bound(left.begin(), left.end(), right[static_cast<std::size_t>(i)]) - left.begin());
    }

    std::vector<int> by_right(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        by_right[static_cast<std::size_t>(i)] = i;
    }
    std::sort(by_right.begin(), by_right.end(), [&](const int a, const int b) {
        return right[static_cast<std::size_t>(a)] < right[static_cast<std::size_t>(b)];
    });

    std::vector<i64> inside(static_cast<std::size_t>(n), 0);
    std::vector<i64> ways(static_cast<std::size_t>(n + 1), 1);
    std::vector<i64> delta(static_cast<std::size_t>(n + 1), 0);
    std::vector<unsigned char> active(static_cast<std::size_t>(n), 0);

    for (const int i : by_right) {
        inside[static_cast<std::size_t>(i)] = ways[static_cast<std::size_t>(i + 1)];

        const i64 inc = mul_mod(inside[static_cast<std::size_t>(i)],
                                ways[static_cast<std::size_t>(next[static_cast<std::size_t>(i)])],
                                mod);
        delta[static_cast<std::size_t>(i)] = inc;
        ways[static_cast<std::size_t>(i)] = add_mod(ways[static_cast<std::size_t>(i)], inc, mod);

        for (int p = i - 1; p >= 0; --p) {
            i64 d = delta[static_cast<std::size_t>(p + 1)];
            if (active[static_cast<std::size_t>(p)] != 0 &&
                next[static_cast<std::size_t>(p)] <= i) {
                d = add_mod(d,
                            mul_mod(inside[static_cast<std::size_t>(p)],
                                    delta[static_cast<std::size_t>(next[static_cast<std::size_t>(p)])],
                                    mod),
                            mod);
            }
            delta[static_cast<std::size_t>(p)] = d;
            ways[static_cast<std::size_t>(p)] = add_mod(ways[static_cast<std::size_t>(p)], d, mod);
        }

        active[static_cast<std::size_t>(i)] = 1;
    }

    return ways[0];
}

bool crosses(const Interval& a, const Interval& b) {
    Interval x = a;
    Interval y = b;
    if (y.left < x.left) {
        std::swap(x, y);
    }
    return x.left < y.left && y.left < x.right && x.right < y.right;
}

i64 brute_connectivity(const std::vector<int>& values) {
    const std::vector<Interval> intervals = build_intervals(values);
    const int n = static_cast<int>(intervals.size());
    assert(n <= 20);

    i64 total = 0;
    for (std::uint64_t mask = 0; mask < (1ULL << n); ++mask) {
        bool ok = true;
        for (int i = 0; i < n && ok; ++i) {
            if (((mask >> i) & 1ULL) == 0ULL) {
                continue;
            }
            for (int j = i + 1; j < n; ++j) {
                if (((mask >> j) & 1ULL) != 0ULL &&
                    crosses(intervals[static_cast<std::size_t>(i)],
                            intervals[static_cast<std::size_t>(j)])) {
                    ok = false;
                    break;
                }
            }
        }
        if (ok) {
            ++total;
        }
    }
    return total;
}

void run_checkpoints() {
    const std::vector<std::vector<int>> cases = {
        {0, 1, 0, 1},
        {0, 0, 1, 2, 2, 1},
        {0, 1, 2, 1, 0, 2},
        {0, 1, 2, 2, 1, 0},
        {0, 1, 2, 3, 1, 4, 0, 5, 4, 2, 6, 7, 3, 8, 6, 5, 9, 8, 9, 7},
    };

    const std::vector<i64> expected = {3, 8, 5, 8, 86};
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const i64 brute = brute_connectivity(cases[i]);
        assert(brute == expected[i]);
        assert(connectivity_number(cases[i], MOD) == brute % MOD);
    }
}

}  // namespace

int main() {
    run_checkpoints();

    const std::vector<int> input = load_csv("resources/documents/1001_input.txt");
    assert(input.size() == 40'000U);

    std::cout << connectivity_number(input, MOD) << '\n';
    return 0;
}
