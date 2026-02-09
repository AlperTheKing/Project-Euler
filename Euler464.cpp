#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

struct Options {
    int n = 20'000'000;
    bool run_checkpoints = true;
};

struct Point {
    int x = 0;
    int y = 0;
    int w = 0;
};

class Fenwick {
public:
    explicit Fenwick(int n) : tree_(static_cast<std::size_t>(n + 1), 0LL) {}

    void add(int idx, const i64 delta) {
        const int n = static_cast<int>(tree_.size()) - 1;
        while (idx <= n) {
            tree_[static_cast<std::size_t>(idx)] += delta;
            idx += idx & -idx;
        }
    }

    i64 sum(int idx) const {
        i64 result = 0;
        while (idx > 0) {
            result += tree_[static_cast<std::size_t>(idx)];
            idx -= idx & -idx;
        }
        return result;
    }

private:
    std::vector<i64> tree_;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n <= 0) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

std::vector<std::int8_t> mobius_values(const int n) {
    std::vector<int> lp(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    std::vector<std::int8_t> mu(static_cast<std::size_t>(n + 1), 0);
    mu[1] = 1;

    for (int i = 2; i <= n; ++i) {
        if (lp[static_cast<std::size_t>(i)] == 0) {
            lp[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (const int p : primes) {
            const int64_t x = static_cast<int64_t>(p) * static_cast<int64_t>(i);
            if (x > n || p > lp[static_cast<std::size_t>(i)]) {
                break;
            }
            lp[static_cast<std::size_t>(x)] = p;
            if (p == lp[static_cast<std::size_t>(i)]) {
                mu[static_cast<std::size_t>(x)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(x)] =
                static_cast<std::int8_t>(-mu[static_cast<std::size_t>(i)]);
        }
    }
    return mu;
}

std::vector<Point> build_compressed_points(const std::vector<std::int8_t>& mu, int& y_shift, int& y_max) {
    std::vector<Point> points;
    points.reserve(mu.size());
    points.push_back(Point{0, 0, 1});  // prefix 0

    int q = 0;
    int m = 0;
    int min_y = 0;
    int max_y = 0;

    for (int i = 1; i < static_cast<int>(mu.size()); ++i) {
        const int v = static_cast<int>(mu[static_cast<std::size_t>(i)]);
        if (v == 0) {
            ++points.back().w;
            continue;
        }

        ++q;
        m += v;
        const int x = q - 199 * m;
        const int y = q + 199 * m;

        if (y < min_y) {
            min_y = y;
        }
        if (y > max_y) {
            max_y = y;
        }

        points.push_back(Point{x, y, 1});
    }

    y_shift = 1 - min_y;
    y_max = max_y + y_shift + 2;
    for (Point& p : points) {
        p.y += y_shift;
    }
    return points;
}

class CDQCounter {
public:
    CDQCounter(std::vector<Point> points, const int fenwick_size)
        : points_(std::move(points)),
          temp_(points_.size()),
          bit_(fenwick_size) {}

    i64 count_pairs() {
        i64 result = 0;
        for (const Point& p : points_) {
            result += static_cast<i64>(p.w) * static_cast<i64>(p.w - 1) / 2LL;
        }
        result_ = result;
        cdq(0, static_cast<int>(points_.size()) - 1);
        return result_;
    }

private:
    std::vector<Point> points_;
    std::vector<Point> temp_;
    Fenwick bit_;
    i64 result_ = 0;

    void cdq(const int left, const int right) {
        if (left >= right) {
            return;
        }
        const int mid = left + (right - left) / 2;
        cdq(left, mid);
        cdq(mid + 1, right);

        int i = left;
        for (int j = mid + 1; j <= right; ++j) {
            while (i <= mid && points_[static_cast<std::size_t>(i)].x <=
                                  points_[static_cast<std::size_t>(j)].x) {
                const Point& p = points_[static_cast<std::size_t>(i)];
                bit_.add(p.y, p.w);
                ++i;
            }
            const Point& r = points_[static_cast<std::size_t>(j)];
            result_ += static_cast<i64>(r.w) * bit_.sum(r.y);
        }
        for (int k = left; k < i; ++k) {
            const Point& p = points_[static_cast<std::size_t>(k)];
            bit_.add(p.y, -static_cast<i64>(p.w));
        }

        int a = left;
        int b = mid + 1;
        int t = left;
        while (a <= mid && b <= right) {
            if (points_[static_cast<std::size_t>(a)].x <= points_[static_cast<std::size_t>(b)].x) {
                temp_[static_cast<std::size_t>(t++)] = points_[static_cast<std::size_t>(a++)];
            } else {
                temp_[static_cast<std::size_t>(t++)] = points_[static_cast<std::size_t>(b++)];
            }
        }
        while (a <= mid) {
            temp_[static_cast<std::size_t>(t++)] = points_[static_cast<std::size_t>(a++)];
        }
        while (b <= right) {
            temp_[static_cast<std::size_t>(t++)] = points_[static_cast<std::size_t>(b++)];
        }
        for (int k = left; k <= right; ++k) {
            points_[static_cast<std::size_t>(k)] = temp_[static_cast<std::size_t>(k)];
        }
    }
};

i64 solve(const int n) {
    const std::vector<std::int8_t> mu = mobius_values(n);

    int y_shift = 0;
    int y_max = 0;
    std::vector<Point> points = build_compressed_points(mu, y_shift, y_max);

    CDQCounter counter(std::move(points), y_max);
    return counter.count_pairs();
}

bool run_checkpoints() {
    if (solve(10) != 13LL) {
        std::cerr << "Checkpoint failed: C(10)\n";
        return false;
    }
    if (solve(500) != 16'676LL) {
        std::cerr << "Checkpoint failed: C(500)\n";
        return false;
    }
    if (solve(10'000) != 20'155'319LL) {
        std::cerr << "Checkpoint failed: C(10000)\n";
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
