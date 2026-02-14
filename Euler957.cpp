#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <algorithm>
#include <functional>

namespace {

using i32 = std::int32_t;
using u128 = unsigned __int128;

u128 pow2_u128(int e) {
    return (u128(1) << e);
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int d = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::pair<u128, u128> day_k_q(int day) {
    if (day == 0) {
        return {2, 2};
    }

    if ((day & 1) == 0) {
        const u128 t = (pow2_u128(day) - 1) / 3;
        const u128 t2 = t * t;
        const u128 t3 = t2 * t;
        const u128 t4 = t2 * t2;
        const u128 k = 3 * t2 + 5 * t + 2;
        const u128 q = (21 * t4 + 70 * t3 + 87 * t2 + 46 * t + 8) / 4;
        return {k, q};
    }

    const u128 t = (pow2_u128(day) - 2) / 3;
    const u128 t2 = t * t;
    const u128 t3 = t2 * t;
    const u128 t4 = t2 * t2;
    const u128 k = 3 * t2 + 7 * t + 4;
    const u128 q = (21 * t4 + 98 * t3 + 171 * t2 + 134 * t + 40) / 4;
    return {k, q};
}

u128 g_formula(int n) {
    if (n == 0) {
        return 2;
    }
    const auto [k, q] = day_k_q(n - 1);
    return 3 * k * k - 2 * q;
}

struct Vec2 {
    i32 x{0};
    i32 y{0};
    bool operator<(const Vec2& other) const {
        if (x != other.x) {
            return x < other.x;
        }
        return y < other.y;
    }
};

Vec2 add(Vec2 a, Vec2 b) {
    return Vec2{static_cast<i32>(a.x + b.x), static_cast<i32>(a.y + b.y)};
}

Vec2 sub(Vec2 a, Vec2 b) {
    return Vec2{static_cast<i32>(a.x - b.x), static_cast<i32>(a.y - b.y)};
}

struct PointState {
    Vec2 a;
    Vec2 b;
    bool operator<(const PointState& other) const {
        if (a.x != other.a.x) {
            return a.x < other.a.x;
        }
        if (a.y != other.a.y) {
            return a.y < other.a.y;
        }
        if (b.x != other.b.x) {
            return b.x < other.b.x;
        }
        return b.y < other.b.y;
    }
};

u128 g_bruteforce_small(int n) {
    std::set<PointState> points;
    points.insert(PointState{Vec2{1, 0}, Vec2{0, 1}});
    points.insert(PointState{Vec2{0, 0}, Vec2{0, 0}});

    for (int day = 0; day < n; ++day) {
        std::set<Vec2> A;
        std::set<Vec2> B;
        std::set<Vec2> C;
        for (const auto& p : points) {
            A.insert(p.a);
            B.insert(p.b);
            C.insert(sub(p.b, p.a));
        }

        std::set<PointState> next;

        for (const auto& a : A) {
            for (const auto& b : B) {
                next.insert(PointState{a, b});
            }
        }
        for (const auto& b : B) {
            for (const auto& c : C) {
                next.insert(PointState{sub(b, c), b});
            }
        }
        for (const auto& a : A) {
            for (const auto& c : C) {
                next.insert(PointState{a, add(a, c)});
            }
        }

        points.swap(next);
    }

    return static_cast<u128>(points.size());
}

void run_validations() {
    assert(g_formula(1) == 8);
    assert(g_formula(2) == 28);
    for (int n = 1; n <= 5; ++n) {
        assert(g_formula(n) == g_bruteforce_small(n));
    }
}

}  // namespace

int main() {
    run_validations();
    std::cout << to_string_u128(g_formula(16)) << '\n';
    return 0;
}
