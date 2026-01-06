#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using int64 = long long;
using ld = long double;

static inline ld clamp_ld(ld x, ld lo, ld hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

static inline int64 safe_floor(ld x) {
    constexpr ld EPS = 1e-13L;
    ld y = x + EPS;
    int64 f = (int64)std::floor(y);

    while ((ld)(f + 1) <= y) ++f;
    while ((ld)f > y) --f;
    return f;
}

static inline int64 g_value(int s, int p, int q, ld PI, ld TWO_PI) {
    const int A = s + p;
    const int B = s + q;

    const ld D = (ld)p + (ld)q - TWO_PI;

    const ld x_beta =
        clamp_ld(((ld)B * B + D * D - (ld)A * A) / (2.0L * (ld)B * D), -1.0L, 1.0L);
    const ld x_delta =
        clamp_ld(((ld)A * A + D * D - (ld)B * B) / (2.0L * (ld)A * D), -1.0L, 1.0L);

    const ld beta = std::acos(x_beta);
    const ld delta = std::acos(x_delta);

    const ld m_max = ((ld)B * beta - (ld)A * delta) / PI;
    const int64 t = safe_floor(m_max);

    const int64 g = (int64)A + t;
    return (g > 0) ? g : 0;
}

static int64 G_value(int n, int threads) {
    const ld PI = std::acos((ld)-1.0L);
    const ld TWO_PI = 2.0L * PI;

    std::atomic<int> next_s{5};
    std::vector<int64> partial((size_t)threads, 0);
    std::vector<std::thread> pool;
    pool.reserve((size_t)threads);

    for (int tid = 0; tid < threads; ++tid) {
        pool.emplace_back([&, tid]() {
            int64 local = 0;
            while (true) {
                int s = next_s.fetch_add(1, std::memory_order_relaxed);
                if (s > n) break;

                const int m = n - s;
                const int max_p = (m - 1) / 2;
                if (max_p < 5) continue;

                for (int p = 5; p <= max_p; ++p) {
                    const int max_q = m - p;
                    for (int q = p + 1; q <= max_q; ++q) {
                        local += g_value(s, p, q, PI, TWO_PI);
                    }
                }
            }
            partial[(size_t)tid] = local;
        });
    }

    for (auto &t : pool) t.join();

    int64 total = 0;
    for (auto v : partial) total += v;
    return total;
}

int main() {
    int threads = (int)std::thread::hardware_concurrency();
    if (threads <= 0) threads = 1;

    assert(G_value(16, 1) == 9);
    assert(G_value(20, 1) == 205);

    {
        const ld PI = std::acos((ld)-1.0L);
        const ld TWO_PI = 2.0L * PI;
        assert(g_value(5, 5, 6, PI, TWO_PI) == 9);
    }

    assert(G_value(100, 1) == 1790942);

    std::cout << G_value(500, threads) << "\n";
    return 0;
}
