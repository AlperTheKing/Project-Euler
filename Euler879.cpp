#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Counter {
    int w;
    int h;
    int n;
    std::vector<int> xs;
    std::vector<int> ys;
    std::vector<u32> between;
    std::vector<u64> memo;

    Counter(int width, int height) : w(width), h(height), n(width * height), xs(n), ys(n) {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                int id = y * w + x;
                xs[id] = x;
                ys[id] = y;
            }
        }

        between.assign(static_cast<std::size_t>(n * n), 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                int dx = xs[j] - xs[i];
                int dy = ys[j] - ys[i];
                int g = std::gcd(std::abs(dx), std::abs(dy));
                if (g <= 1) continue;
                int stepx = dx / g;
                int stepy = dy / g;
                u32 mask = 0;
                int cx = xs[i];
                int cy = ys[i];
                for (int t = 1; t < g; ++t) {
                    cx += stepx;
                    cy += stepy;
                    int mid = cy * w + cx;
                    mask |= (1U << mid);
                }
                between[static_cast<std::size_t>(i * n + j)] = mask;
            }
        }

        std::size_t total_states = static_cast<std::size_t>(1U << n) * static_cast<std::size_t>(n);
        memo.assign(total_states, std::numeric_limits<u64>::max());
    }

    u64 dfs(u32 used, int last) {
        std::size_t idx = static_cast<std::size_t>(used) * static_cast<std::size_t>(n) + static_cast<std::size_t>(last);
        u64& cache = memo[idx];
        if (cache != std::numeric_limits<u64>::max()) return cache;

        u64 ans = 0;
        u32 free_mask = ((1U << n) - 1U) ^ used;

        for (int nxt = 0; nxt < n; ++nxt) {
            u32 bit = (1U << nxt);
            if ((free_mask & bit) == 0) continue;
            u32 need = between[static_cast<std::size_t>(last * n + nxt)];
            if ((need & free_mask) != 0) continue;
            u32 used2 = used | bit;
            ans += 1 + dfs(used2, nxt);
        }

        cache = ans;
        return ans;
    }

    u64 count_passwords() {
        u64 total = 0;
        for (int start = 0; start < n; ++start) {
            total += dfs(1U << start, start);
        }
        return total;
    }
};

int main() {
    {
        Counter c3(3, 3);
        assert(c3.count_passwords() == 389'488ULL);
    }

    Counter c4(4, 4);
    std::cout << c4.count_passwords() << '\n';
    return 0;
}
