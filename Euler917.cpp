#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 kMod = 998388889LL;
constexpr i64 kSeed = 102022661LL;

struct ABData {
    std::vector<int> a;
    std::vector<int> b;
    i64 sum_a;
    i64 sum_b;
};

ABData generate_ab(int n) {
    ABData data;
    data.a.assign(n + 1, 0);
    data.b.assign(n + 1, 0);
    data.sum_a = 0;
    data.sum_b = 0;

    i64 s = kSeed;
    for (int i = 1; i <= n; ++i) {
        data.a[i] = static_cast<int>(s);
        data.sum_a += s;

        s = static_cast<i64>((__int128)s * s % kMod);
        data.b[i] = static_cast<int>(s);
        data.sum_b += s;

        if (i != n) {
            s = static_cast<i64>((__int128)s * s % kMod);
        }
    }
    return data;
}

std::vector<int> lower_hull_indices(const std::vector<int>& b, int max_t) {
    std::vector<int> hull;
    hull.reserve(64);

    auto cross = [&](int t1, int t2, int t3) -> __int128 {
        const __int128 x1 = t2 - t1;
        const __int128 y1 = static_cast<i64>(b[t2 + 1]) - static_cast<i64>(b[t1 + 1]);
        const __int128 x2 = t3 - t1;
        const __int128 y2 = static_cast<i64>(b[t3 + 1]) - static_cast<i64>(b[t1 + 1]);
        return x1 * y2 - y1 * x2;
    };

    for (int t = 0; t <= max_t; ++t) {
        while (hull.size() >= 2U) {
            const int t1 = hull[hull.size() - 2];
            const int t2 = hull[hull.size() - 1];
            if (cross(t1, t2, t) <= 0) {
                hull.pop_back();
            } else {
                break;
            }
        }
        hull.push_back(t);
    }

    return hull;
}

i64 minimal_extra_cost(const std::vector<int>& a, const std::vector<int>& b, int n) {
    const int edge_steps = n - 1;
    if (edge_steps == 0) {
        return 0;
    }

    const std::vector<int> hull_t = lower_hull_indices(b, edge_steps);
    const int m = static_cast<int>(hull_t.size());

    std::vector<i64> t_vals(m);
    std::vector<i64> g_vals(m);
    for (int k = 0; k < m; ++k) {
        t_vals[k] = hull_t[k];
        g_vals[k] = static_cast<i64>(b[hull_t[k] + 1]);
    }

    std::vector<i64> prev(m, 0);
    std::vector<i64> cur(m, 0);

    for (int i = 1; i <= edge_steps; ++i) {
        const i64 u = static_cast<i64>(a[i]) - static_cast<i64>(a[i + 1]);
        i64 pref = std::numeric_limits<i64>::max();
        for (int k = 0; k < m; ++k) {
            if (prev[k] < pref) {
                pref = prev[k];
            }
            cur[k] = pref + u * t_vals[k] + g_vals[k];
        }
        prev.swap(cur);
    }

    i64 best = prev[0];
    for (int k = 1; k < m; ++k) {
        if (prev[k] < best) {
            best = prev[k];
        }
    }

    return static_cast<i64>(edge_steps) * static_cast<i64>(a[n]) + best;
}

i64 solve_fast(int n) {
    const ABData data = generate_ab(n);
    const i64 extra = minimal_extra_cost(data.a, data.b, n);
    return data.sum_a + data.sum_b + extra;
}

i64 solve_quadratic_check(int n) {
    const ABData data = generate_ab(n);
    const i64 inf = std::numeric_limits<i64>::max() / 4;
    std::vector<i64> prev(n + 1, inf);
    std::vector<i64> cur(n + 1, inf);

    prev[1] = 0;
    for (int i = 1; i <= n; ++i) {
        std::fill(cur.begin(), cur.end(), inf);
        for (int j = 1; j <= n; ++j) {
            if (i == 1 && j == 1) {
                cur[j] = 0;
                continue;
            }
            i64 best = inf;
            if (j > 1 && cur[j - 1] != inf) {
                best = std::min(best, cur[j - 1] + static_cast<i64>(data.a[i]));
            }
            if (prev[j] != inf) {
                best = std::min(best, prev[j] + static_cast<i64>(data.b[j]));
            }
            cur[j] = best;
        }
        prev.swap(cur);
    }

    return data.sum_a + data.sum_b + prev[n];
}

void validate() {
    assert(solve_fast(1) == 966774091LL);
    assert(solve_fast(2) == 2388327490LL);
    assert(solve_fast(10) == 13389278727LL);
    for (int n = 1; n <= 120; ++n) {
        assert(solve_fast(n) == solve_quadratic_check(n));
    }
}

}  // namespace

int main() {
    validate();
    std::cout << solve_fast(10'000'000) << '\n';
    return 0;
}
