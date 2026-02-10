#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;

static constexpr int C_BITS = 23;
static constexpr int B_BITS = 12;
static constexpr int A_BITS = 12;
static constexpr i64 C_OFF = 4'000'000;

static inline u64 pack_line(int A, int B, i64 C) {
    const i64 Co = C + C_OFF;
    assert(0 <= A && A < (1 << A_BITS));
    assert(-2048 <= B && B < 2048);
    assert(0 <= Co && Co < (1LL << C_BITS));
    const u64 Ao = (u64)A;
    const u64 Bo = (u64)(B + 2048);
    const u64 Co_u = (u64)Co;
    return (Ao << (B_BITS + C_BITS)) | (Bo << C_BITS) | Co_u;
}

static std::vector<std::pair<int, int>> generate_points(int n) {
    std::vector<std::pair<int, int>> pts;
    pts.reserve(n);
    i64 S = 290797;
    for (int k = 0; k < n; ++k) {
        S = (S * S) % 50515093;
        const int x = (int)(S % 2000) - 1000;
        S = (S * S) % 50515093;
        const int y = (int)(S % 2000) - 1000;
        pts.emplace_back(x, y);
    }
    return pts;
}

static std::pair<i64, i64> compute_MS(int n) {
    const auto pts = generate_points(n);
    std::vector<u64> keys;
    keys.reserve((i64)n * (n - 1) / 2);

    for (int i = 0; i < n; ++i) {
        const auto [x1, y1] = pts[i];
        for (int j = i + 1; j < n; ++j) {
            const auto [x2, y2] = pts[j];
            int dx = x2 - x1;
            int dy = y2 - y1;
            if (dx == 0 && dy == 0) continue;
            const int g = std::gcd(std::abs(dx), std::abs(dy));
            dx /= g;
            dy /= g;

            int A = dy;
            int B = -dx;
            i64 C = (i64)A * x1 + (i64)B * y1;
            if (A < 0 || (A == 0 && B < 0)) {
                A = -A;
                B = -B;
                C = -C;
            }
            keys.push_back(pack_line(A, B, C));
        }
    }

    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    const i64 M = (i64)keys.size();
    i64 parallel2 = 0;
    for (i64 i = 0; i < M;) {
        const u64 slope = keys[(size_t)i] >> C_BITS;
        i64 j = i + 1;
        while (j < M && (keys[(size_t)j] >> C_BITS) == slope) ++j;
        const i64 m = j - i;
        parallel2 += m * (m - 1);
        i = j;
    }
    const i64 S = M * (M - 1) - parallel2;
    return {M, S};
}

int main() {
    {
        const auto pts3 = generate_points(3);
        assert((pts3[0] == std::make_pair(527, 144)));
        assert((pts3[1] == std::make_pair(-488, 732)));
        assert((pts3[2] == std::make_pair(-454, -947)));
    }

    {
        const auto [m3, s3] = compute_MS(3);
        assert(m3 == 3);
        assert(s3 == 6);
    }

    {
        const auto [m100, s100] = compute_MS(100);
        assert(m100 == 4948);
        assert(s100 == 24477690);
    }

    std::cout << compute_MS(2500).second << "\n";
    return 0;
}

