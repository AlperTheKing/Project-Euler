#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr int TARGET = 1'000'000;

struct Leg {
    int offset;
    int hypotenuse;
};

i64 isqrt(const i64 n) {
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) {
        ++r;
    }
    while (r * r > n) {
        --r;
    }
    return r;
}

bool is_square(const i64 n, i64& root) {
    root = isqrt(n);
    return root * root == n;
}

long double normalize_angle(long double angle) {
    const long double half_pi = std::acos(-1.0L) / 2.0L;
    angle = std::fmod(angle, half_pi);
    if (angle < 0) {
        angle += half_pi;
    }
    if (angle < 1e-18L || half_pi - angle < 1e-18L) {
        return 0.0L;
    }
    return angle;
}

struct WidthState {
    long double width;
    int u_min;
    int u_max;
    int v_min;
    int v_max;
};

long double minimum_square(int a, int b, int c) {
    if (a > b) {
        std::swap(a, b);
    }
    if (b > c) {
        std::swap(b, c);
    }
    if (a > b) {
        std::swap(a, b);
    }

    const long double x = (static_cast<long double>(a) * a + static_cast<long double>(c) * c -
                           static_cast<long double>(b) * b) /
                          (2.0L * c);
    long double y2 = static_cast<long double>(a) * a - x * x;
    if (y2 < 0 && y2 > -1e-12L) {
        y2 = 0;
    }
    const long double y = std::sqrt(y2);
    const std::array<std::pair<long double, long double>, 3> points{{{0, 0}, {static_cast<long double>(c), 0}, {x, y}}};

    auto evaluate = [&](const long double angle) {
        const long double co = std::cos(angle);
        const long double si = std::sin(angle);
        std::array<long double, 3> u{};
        std::array<long double, 3> v{};
        for (int i = 0; i < 3; ++i) {
            u[static_cast<std::size_t>(i)] = points[static_cast<std::size_t>(i)].first * co +
                                             points[static_cast<std::size_t>(i)].second * si;
            v[static_cast<std::size_t>(i)] = -points[static_cast<std::size_t>(i)].first * si +
                                             points[static_cast<std::size_t>(i)].second * co;
        }

        const int u_min = static_cast<int>(std::min_element(u.begin(), u.end()) - u.begin());
        const int u_max = static_cast<int>(std::max_element(u.begin(), u.end()) - u.begin());
        const int v_min = static_cast<int>(std::min_element(v.begin(), v.end()) - v.begin());
        const int v_max = static_cast<int>(std::max_element(v.begin(), v.end()) - v.begin());
        return WidthState{std::max(u[static_cast<std::size_t>(u_max)] - u[static_cast<std::size_t>(u_min)],
                                   v[static_cast<std::size_t>(v_max)] - v[static_cast<std::size_t>(v_min)]),
                          u_min,
                          u_max,
                          v_min,
                          v_max};
    };

    const long double half_pi = std::acos(-1.0L) / 2.0L;
    std::vector<long double> angles{0};
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const long double edge = std::atan2(points[static_cast<std::size_t>(j)].second -
                                                    points[static_cast<std::size_t>(i)].second,
                                                points[static_cast<std::size_t>(j)].first -
                                                    points[static_cast<std::size_t>(i)].first);
            angles.push_back(normalize_angle(edge));
            angles.push_back(normalize_angle(edge + half_pi));
        }
    }

    std::sort(angles.begin(), angles.end());
    std::vector<long double> boundaries;
    for (const long double angle : angles) {
        if (boundaries.empty() || std::fabs(angle - boundaries.back()) > 1e-15L) {
            boundaries.push_back(angle);
        }
    }

    std::vector<long double> candidates = boundaries;
    std::vector<long double> endpoints = boundaries;
    endpoints.push_back(half_pi);
    for (std::size_t i = 0; i + 1 < endpoints.size(); ++i) {
        const long double lo = endpoints[i];
        const long double hi = endpoints[i + 1];
        if (hi - lo < 1e-15L) {
            continue;
        }

        const WidthState state = evaluate((lo + hi) / 2.0L);
        const auto& a0 = points[static_cast<std::size_t>(state.u_min)];
        const auto& a1 = points[static_cast<std::size_t>(state.u_max)];
        const auto& b0 = points[static_cast<std::size_t>(state.v_min)];
        const auto& b1 = points[static_cast<std::size_t>(state.v_max)];
        const long double ax = a1.first - a0.first;
        const long double ay = a1.second - a0.second;
        const long double bx = b1.first - b0.first;
        const long double by = b1.second - b0.second;
        const long double left = ax - by;
        const long double right = ay + bx;
        if (std::fabs(right) < 1e-18L) {
            continue;
        }

        const long double root = std::atan2(-left, right);
        for (const long double angle : {root, root + std::acos(-1.0L)}) {
            const long double current = normalize_angle(angle);
            if (lo + 1e-15L < current && current < hi - 1e-15L) {
                candidates.push_back(current);
            }
        }
    }

    long double best = std::numeric_limits<long double>::max();
    for (const long double angle : candidates) {
        best = std::min(best, evaluate(angle).width);
    }
    return best;
}

u64 triangle_key(int a, int b, int c) {
    if (a > b) {
        std::swap(a, b);
    }
    if (b > c) {
        std::swap(b, c);
    }
    if (a > b) {
        std::swap(a, b);
    }
    return (static_cast<u64>(a) << 42) | (static_cast<u64>(b) << 21) | static_cast<u64>(c);
}

std::vector<std::vector<Leg>> pythagorean_offsets(const int limit) {
    std::vector<std::vector<Leg>> by_side(static_cast<std::size_t>(limit + 1));
    const int bound = isqrt(2LL * limit) + 2;

    for (int m = 2; m <= bound; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0 || std::gcd(m, n) != 1) {
                continue;
            }

            const int a = m * m - n * n;
            const int b = 2 * m * n;
            const int c = m * m + n * n;

            for (const auto [side0, offset0] : {std::pair<int, int>{a, b}, {b, a}}) {
                for (int k = 1; k * side0 <= limit; ++k) {
                    const int side = k * side0;
                    const int offset = k * offset0;
                    if (offset <= side) {
                        by_side[static_cast<std::size_t>(side)].push_back({offset, k * c});
                    }
                }
            }
        }
    }

    for (std::vector<Leg>& legs : by_side) {
        std::sort(legs.begin(), legs.end(), [](const Leg& lhs, const Leg& rhs) {
            if (lhs.offset != rhs.offset) {
                return lhs.offset < rhs.offset;
            }
            return lhs.hypotenuse < rhs.hypotenuse;
        });
        legs.erase(std::unique(legs.begin(), legs.end(), [](const Leg& lhs, const Leg& rhs) {
                       return lhs.offset == rhs.offset && lhs.hypotenuse == rhs.hypotenuse;
                   }),
                   legs.end());
    }

    return by_side;
}

u64 solve(const int limit) {
    const std::vector<std::vector<Leg>> legs = pythagorean_offsets(limit);
    std::unordered_set<u64> seen;
    u64 total = 0;

    auto add_triangle = [&](const int a, const int b, const int c) {
        const u64 key = triangle_key(a, b, c);
        if (seen.insert(key).second) {
            total += static_cast<u64>(a + b + c);
        }
    };

    for (int side = 1; side <= limit; ++side) {
        const std::vector<Leg>& current = legs[static_cast<std::size_t>(side)];
        for (std::size_t i = 0; i < current.size(); ++i) {
            const i64 p = current[i].offset;
            for (std::size_t j = i; j < current.size(); ++j) {
                const i64 q = current[j].offset;
                const i64 base = p + q;
                if (base <= side && p * q >= static_cast<i64>(side) * (side - base)) {
                    add_triangle(static_cast<int>(base), current[i].hypotenuse, current[j].hypotenuse);
                }

                const i64 u = side - p;
                const i64 v = side - q;
                i64 third = 0;
                if (u > 0 && v > 0 && is_square(u * u + v * v, third) &&
                    minimum_square(current[i].hypotenuse, current[j].hypotenuse, static_cast<int>(third)) + 1e-6L >=
                        side) {
                    add_triangle(current[i].hypotenuse, current[j].hypotenuse, static_cast<int>(third));
                }
            }
        }
    }

    return total;
}

void run_checkpoints() {
    assert(solve(40) == 346);
    assert(solve(400) == 76'402);
    assert(solve(2'000) == 3'237'036);
}

}  // namespace

int main() {
    run_checkpoints();
    std::cout << solve(TARGET) << '\n';
    return 0;
}
