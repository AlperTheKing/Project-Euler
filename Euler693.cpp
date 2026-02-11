#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <queue>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

u32 l_xy(const u32 x, const u32 y) {
    u32 z = x;
    u32 a = y;
    u32 length = 1;

    while (a > 1U) {
        a = static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(a)) % z);
        ++z;
        ++length;
    }

    return length;
}

std::vector<u32> initial_active_after_first_step(const u32 x) {
    if (x <= 2U) {
        return {};
    }

    std::vector<std::uint8_t> seen(x, 0U);
    std::vector<u32> active;
    active.reserve(x / 4U);

    const u64 mod = static_cast<u64>(x);
    const u64 limit = mod / 2ULL;

    u64 y = 2ULL;
    u64 sq = (y * y) % mod;
    u64 delta = 2ULL * y + 1ULL;

    while (y <= limit) {
        if (sq > 1ULL && !seen[static_cast<std::size_t>(sq)]) {
            seen[static_cast<std::size_t>(sq)] = 1U;
            active.push_back(static_cast<u32>(sq));
        }

        sq += delta;
        if (sq >= mod) {
            sq -= mod;
            if (sq >= mod) {
                sq -= mod;
            }
        }
        delta += 2ULL;
        ++y;
    }

    return active;
}

std::vector<u32> step_active_dense(const std::vector<u32>& active, const u32 mod) {
    std::vector<std::uint8_t> seen(mod, 0U);
    std::vector<u32> next;
    next.reserve(active.size());

    for (const u32 a : active) {
        const u32 v = static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(a)) % mod);
        if (v > 1U && !seen[static_cast<std::size_t>(v)]) {
            seen[static_cast<std::size_t>(v)] = 1U;
            next.push_back(v);
        }
    }

    return next;
}

std::vector<u32> step_active_sparse(const std::vector<u32>& active, const u32 mod) {
    std::unordered_set<u32> next_set;
    next_set.reserve(active.size() * 2ULL);

    for (const u32 a : active) {
        const u32 v = static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(a)) % mod);
        if (v > 1U) {
            next_set.insert(v);
        }
    }

    std::vector<u32> next;
    next.reserve(next_set.size());
    for (const u32 v : next_set) {
        next.push_back(v);
    }
    return next;
}

u32 g_value(const u32 x, const std::size_t dense_threshold = 120000U) {
    if (x < 2U) {
        return 0U;
    }
    if (x == 2U) {
        return 1U;
    }

    std::vector<u32> active = initial_active_after_first_step(x);

    u32 length = 2U;
    u32 mod = x + 1U;

    while (active.size() > 1U) {
        if (active.size() >= dense_threshold) {
            active = step_active_dense(active, mod);
        } else {
            active = step_active_sparse(active, mod);
        }
        ++length;
        ++mod;
        if (active.empty()) {
            return length;
        }
    }

    if (active.empty()) {
        return length;
    }

    u32 value = active[0];
    while (value > 1U) {
        value = static_cast<u32>((static_cast<u64>(value) * static_cast<u64>(value)) % mod);
        ++length;
        ++mod;
    }

    return length;
}

u32 snapped_grid_step(const u32 raw_step) {
    if (raw_step <= 1U) {
        return 1U;
    }

    const int power = static_cast<int>(std::floor(std::log10(static_cast<double>(raw_step))));
    u32 base = 1U;
    for (int i = 0; i < power; ++i) {
        base *= 10U;
    }

    for (const u32 m : {1U, 2U, 5U, 10U}) {
        if (m * base >= raw_step) {
            return m * base;
        }
    }
    return raw_step;
}

struct Interval {
    u32 upper_bound;
    u32 left;
    u32 right;
    u32 g_right;
};

struct IntervalOrder {
    bool operator()(const Interval& a, const Interval& b) const {
        return a.upper_bound < b.upper_bound;
    }
};

u32 f_value(const u32 n, const u32 target_points = 16U) {
    if (n < 2U) {
        return 0U;
    }

    const u32 points_count = std::max(2U, target_points);
    const u32 raw_step = std::max(1U, (n - 2U) / (points_count - 1U));
    const u32 grid_step = snapped_grid_step(raw_step);

    std::vector<u32> points;
    points.reserve(points_count + 2U);
    for (u32 x = 2U; x <= n; x += grid_step) {
        points.push_back(x);
        if (x > n - grid_step) {
            break;
        }
    }
    if (points.back() != n) {
        points.push_back(n);
    }

    std::vector<u32> g_points(points.size(), 0U);
    {
        const unsigned hw = std::thread::hardware_concurrency();
        const unsigned thread_count =
            std::max(1U, std::min<unsigned>(hw == 0U ? 1U : hw, points.size()));

        std::atomic<std::size_t> next_index(0U);
        std::vector<std::thread> pool;
        pool.reserve(thread_count);

        for (unsigned t = 0U; t < thread_count; ++t) {
            pool.emplace_back([&]() {
                while (true) {
                    const std::size_t i = next_index.fetch_add(1U, std::memory_order_relaxed);
                    if (i >= points.size()) {
                        break;
                    }
                    g_points[i] = g_value(points[i]);
                }
            });
        }

        for (auto& th : pool) {
            th.join();
        }
    }

    std::unordered_map<u32, u32> cache;
    cache.reserve(points.size() * 4ULL);

    u32 best = 0U;
    for (std::size_t i = 0; i < points.size(); ++i) {
        cache[points[i]] = g_points[i];
        best = std::max(best, g_points[i]);
    }

    auto G = [&](const u32 x) -> u32 {
        const auto it = cache.find(x);
        if (it != cache.end()) {
            return it->second;
        }
        const u32 gx = g_value(x);
        cache.emplace(x, gx);
        return gx;
    };

    std::priority_queue<Interval, std::vector<Interval>, IntervalOrder> pq;
    for (std::size_t i = 1; i < points.size(); ++i) {
        const u32 left = points[i - 1U];
        const u32 right = points[i];
        const u32 g_right = g_points[i];
        const u32 upper_bound = g_right + (right - left);
        pq.push(Interval{upper_bound, left, right, g_right});
    }

    while (!pq.empty()) {
        const Interval top = pq.top();
        if (top.upper_bound <= best) {
            break;
        }
        pq.pop();

        if (top.right - top.left <= 1U) {
            continue;
        }

        const u32 mid = (top.left + top.right) / 2U;
        const u32 g_mid = G(mid);
        best = std::max(best, g_mid);

        if (mid - top.left > 1U) {
            const u32 ub_left = g_mid + (mid - top.left);
            if (ub_left > best) {
                pq.push(Interval{ub_left, top.left, mid, g_mid});
            }
        }

        if (top.right - mid > 1U) {
            const u32 ub_right = top.g_right + (top.right - mid);
            if (ub_right > best) {
                pq.push(Interval{ub_right, mid, top.right, top.g_right});
            }
        }
    }

    return best;
}

}  // namespace

int main() {
    assert(l_xy(5U, 3U) == 29U);
    assert(g_value(5U) == 29U);
    assert(f_value(100U) == 145U);
    assert(f_value(10'000U) == 8'824U);

    std::cout << f_value(3'000'000U) << '\n';
    return 0;
}

