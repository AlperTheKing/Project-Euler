#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kMod = 1'000'000'007U;

struct Steps {
    std::vector<int> vertical;
    std::vector<std::pair<int, int>> other;
    int max_dx = 0;
};

std::vector<int> fibonacci_upto(int limit) {
    std::vector<int> fib{1, 2};
    while (true) {
        const int nxt = fib[fib.size() - 1] + fib[fib.size() - 2];
        if (nxt > limit) break;
        fib.push_back(nxt);
    }
    return fib;
}

Steps build_steps(int w, int h) {
    const int max_len = static_cast<int>(std::sqrt(static_cast<long double>(w) * w +
                                                   static_cast<long double>(h) * h));
    const auto fib = fibonacci_upto(max_len);

    std::vector<std::pair<int, int>> all;
    for (int f : fib) {
        const int fsq = f * f;
        const int x_max = std::min(f, w);
        for (int x = 0; x <= x_max; ++x) {
            const int y2 = fsq - x * x;
            if (y2 < 0) continue;
            const int y = static_cast<int>(std::sqrt(static_cast<long double>(y2)));
            if (y * y != y2) continue;
            if (y > h) continue;
            if (x == 0 && y == 0) continue;
            all.emplace_back(x, y);
        }
    }

    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());

    Steps s;
    for (const auto& [dx, dy] : all) {
        if (dx == 0) {
            s.vertical.push_back(dy);
        } else {
            s.other.emplace_back(dx, dy);
            if (dx > s.max_dx) s.max_dx = dx;
        }
    }

    std::sort(s.vertical.begin(), s.vertical.end());
    std::sort(s.other.begin(), s.other.end());
    return s;
}

u32 solve_case(int w, int h) {
    const Steps steps = build_steps(w, h);
    const int row_len = h + 1;
    const int ring = steps.max_dx + 1;

    std::vector<u32> rows(static_cast<std::size_t>(ring) * static_cast<std::size_t>(row_len), 0);
    std::vector<u32> acc(static_cast<std::size_t>(row_len), 0U);

    auto row_ptr = [&](int ridx) -> u32* {
        return rows.data() + static_cast<std::size_t>(ridx) * static_cast<std::size_t>(row_len);
    };

    for (int x = 0; x <= w; ++x) {
        std::fill(acc.begin(), acc.end(), 0U);
        u32* cur = row_ptr(x % ring);
        std::fill(cur, cur + row_len, 0U);

        for (const auto& [dx, dy] : steps.other) {
            if (dx > x) break;
            const u32* src = row_ptr((x - dx) % ring);
            for (int y = dy; y <= h; ++y) {
                u32 v = acc[static_cast<std::size_t>(y)] + src[static_cast<std::size_t>(y - dy)];
                if (v >= kMod) {
                    v -= kMod;
                }
                acc[static_cast<std::size_t>(y)] = v;
            }
        }

        for (int y = 0; y <= h; ++y) {
            u32 val = acc[static_cast<std::size_t>(y)];
            if (x == 0 && y == 0) {
                val += 1;
                if (val >= kMod) val -= kMod;
            }
            for (int dy : steps.vertical) {
                if (dy > y) break;
                val += cur[static_cast<std::size_t>(y - dy)];
                if (val >= kMod) val -= kMod;
            }
            cur[static_cast<std::size_t>(y)] = static_cast<u32>(val);
        }
    }

    const u32* target_row = row_ptr(w % ring);
    return target_row[h];
}

}  // namespace

int main() {
    assert(solve_case(3, 4) == 278U);
    assert(solve_case(10, 10) == 215846462U);
    std::cout << solve_case(10'000, 10'000) << "\n";
    return 0;
}
