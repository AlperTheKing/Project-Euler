#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <vector>
#include <array>
#include <functional>

struct Point {
    int x;
    int y;
};

struct ByY {
    bool operator()(const std::array<int, 3>& a, const std::array<int, 3>& b) const {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2];
    }
};

static std::vector<Point> generate_points(int k) {
    static constexpr std::int64_t mod = 50'515'093;
    std::int64_t s = 290'797;

    std::vector<Point> pts;
    pts.reserve(k);

    for (int i = 0; i < k; ++i) {
        Point p;
        p.x = static_cast<int>(s);
        s = (s * s) % mod;
        p.y = static_cast<int>(s);
        s = (s * s) % mod;
        pts.push_back(p);
    }

    return pts;
}

static inline std::int64_t dist2(const Point& a, const Point& b) {
    const std::int64_t dx = static_cast<std::int64_t>(a.x) - static_cast<std::int64_t>(b.x);
    const std::int64_t dy = static_cast<std::int64_t>(a.y) - static_cast<std::int64_t>(b.y);
    return dx * dx + dy * dy;
}

static long double closest_distance(int k) {
    std::vector<Point> pts = generate_points(k);
    std::sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y;
    });

    std::int64_t best2 = dist2(pts[0], pts[1]);

    std::set<std::array<int, 3>, ByY> active;
    active.insert({pts[0].y, pts[0].x, 0});

    int left = 0;

    for (int i = 1; i < k; ++i) {
        const long double d = std::sqrt(static_cast<long double>(best2));

        while (left < i) {
            const std::int64_t dx = static_cast<std::int64_t>(pts[i].x) - static_cast<std::int64_t>(pts[left].x);
            if (static_cast<long double>(dx) <= d) {
                break;
            }
            active.erase({pts[left].y, pts[left].x, left});
            ++left;
        }

        const int y_lo = static_cast<int>(std::floor(static_cast<long double>(pts[i].y) - d));
        const int y_hi = static_cast<int>(std::ceil(static_cast<long double>(pts[i].y) + d));

        auto it = active.lower_bound({y_lo, std::numeric_limits<int>::min(), std::numeric_limits<int>::min()});
        while (it != active.end() && (*it)[0] <= y_hi) {
            const int idx = (*it)[2];
            best2 = std::min(best2, dist2(pts[i], pts[idx]));
            ++it;
        }

        active.insert({pts[i].y, pts[i].x, i});
    }

    return std::sqrt(static_cast<long double>(best2));
}

int main() {
    const long double d14 = closest_distance(14);
    assert(std::fabsl(d14 - 546446.466846479L) < 1e-9L);

    const long double ans = closest_distance(2'000'000);
    std::cout << std::fixed << std::setprecision(9) << ans << '\n';
    return 0;
}
