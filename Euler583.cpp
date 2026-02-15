#include <pthread.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unistd.h>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Tri {
    int a;
    int t;
    int s;
};

struct Rect {
    int a;
    int h;
};

struct GroupRange {
    int a;
    std::size_t tri_begin;
    std::size_t tri_end;
    std::size_t rect_begin;
    std::size_t rect_end;
};

static std::vector<Tri> generate_triangles(int L) {
    std::vector<Tri> tris;

    const int mmax = (int)std::sqrt((long double)L) + 2;
    for (int m = 2; m <= mmax; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0) continue;
            if (std::gcd(m, n) != 1) continue;

            const int a0 = m * m - n * n;
            const int b0 = 2 * m * n;
            const int c0 = m * m + n * n;
            if (c0 > L) continue;

            const int kmax = L / c0;
            for (int k = 1; k <= kmax; ++k) {
                const int a = k * a0;
                const int t = k * b0;
                const int s = k * c0;
                tris.push_back(Tri{a, t, s});
                tris.push_back(Tri{t, a, s});
            }
        }
    }

    return tris;
}

static std::vector<Rect> generate_rectangles(int L) {
    const int max_leg = 2 * L;
    std::vector<Rect> rects;

    const int mmax = (int)std::sqrt((long double)max_leg) + 2;
    for (int m = 2; m <= mmax; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0) continue;
            if (std::gcd(m, n) != 1) continue;

            const int a0 = m * m - n * n;
            const int b0 = 2 * m * n;
            const int max0 = std::max(a0, b0);
            if (max0 > max_leg) continue;

            const int kmax = max_leg / max0;
            for (int k = 1; k <= kmax; ++k) {
                const int x = k * a0;
                const int y = k * b0;

                if ((x & 1) == 0) {
                    const int a = x / 2;
                    const int h = y;
                    if (a > 0 && a <= L && h > 0 && h <= L) rects.push_back(Rect{a, h});
                }
                if ((y & 1) == 0) {
                    const int a = y / 2;
                    const int h = x;
                    if (a > 0 && a <= L && h > 0 && h <= L) rects.push_back(Rect{a, h});
                }
            }
        }
    }

    return rects;
}

static u64 compute_group_sum(const GroupRange& g,
                             const std::vector<Tri>& tris,
                             const std::vector<Rect>& rects,
                             const int L,
                             std::vector<u32>& mark,
                             u32& stamp) {
    if (++stamp == 0U) {
        std::fill(mark.begin(), mark.end(), 0U);
        stamp = 1U;
    }
    for (std::size_t t = g.tri_begin; t < g.tri_end; ++t) {
        const int tv = tris[t].t;
        if (tv >= 0 && tv <= L) {
            mark[static_cast<std::size_t>(tv)] = stamp;
        }
    }

    const int a = g.a;
    u64 sum = 0ULL;
    for (std::size_t r = g.rect_begin; r < g.rect_end; ++r) {
        const int h = rects[r].h;
        if (h <= 0) continue;
        const int limit_s = L - a - h;
        if (limit_s <= 0) break;

        for (std::size_t t = g.tri_begin; t < g.tri_end; ++t) {
            const int tv = tris[t].t;
            if (tv >= h) break;
            const int s = tris[t].s;
            if (s > limit_s) break;

            const int u = h + tv;
            if (u > L) continue;
            if (mark[static_cast<std::size_t>(u)] != stamp) continue;

            sum += 2ULL * (u64)(a + h + s);
        }
    }
    return sum;
}

static int detect_thread_count(std::size_t work_items) {
    long cores = ::sysconf(_SC_NPROCESSORS_ONLN);
    int threads = (cores > 0) ? static_cast<int>(cores) : 4;
    if (threads < 1) threads = 1;
    if (work_items == 0) return 1;
    if ((std::size_t)threads > work_items) threads = static_cast<int>(work_items);
    return threads;
}

struct WorkerTask {
    const std::vector<Tri>* tris = nullptr;
    const std::vector<Rect>* rects = nullptr;
    const std::vector<GroupRange>* groups = nullptr;
    int L = 0;
    std::atomic<std::size_t>* next_idx = nullptr;
    std::vector<u32> mark;
    u32 stamp = 0U;
    u64 partial = 0ULL;
};

static void* worker_entry(void* raw) {
    auto* task = static_cast<WorkerTask*>(raw);
    u64 local = 0ULL;
    if (task->mark.empty()) {
        task->mark.resize(static_cast<std::size_t>(task->L + 1), 0U);
        task->stamp = 0U;
    }
    while (true) {
        const std::size_t idx = task->next_idx->fetch_add(1, std::memory_order_relaxed);
        if (idx >= task->groups->size()) break;
        local += compute_group_sum((*task->groups)[idx], *task->tris, *task->rects, task->L,
                                   task->mark, task->stamp);
    }
    task->partial = local;
    return nullptr;
}

static u64 S(u64 p) {
    const int L = (int)(p / 2);
    const auto tris_raw = generate_triangles(L);
    const auto rects_raw = generate_rectangles(L);

    std::vector<u32> tri_count(static_cast<std::size_t>(L + 1), 0U);
    std::vector<u32> rect_count(static_cast<std::size_t>(L + 1), 0U);
    for (const auto& t : tris_raw) {
        if (t.a >= 1 && t.a <= L) {
            ++tri_count[static_cast<std::size_t>(t.a)];
        }
    }
    for (const auto& r : rects_raw) {
        if (r.a >= 1 && r.a <= L) {
            ++rect_count[static_cast<std::size_t>(r.a)];
        }
    }

    std::vector<u32> tri_off(static_cast<std::size_t>(L + 2), 0U);
    std::vector<u32> rect_off(static_cast<std::size_t>(L + 2), 0U);
    for (int a = 1; a <= L; ++a) {
        tri_off[static_cast<std::size_t>(a + 1)] =
            tri_off[static_cast<std::size_t>(a)] + tri_count[static_cast<std::size_t>(a)];
        rect_off[static_cast<std::size_t>(a + 1)] =
            rect_off[static_cast<std::size_t>(a)] + rect_count[static_cast<std::size_t>(a)];
    }

    std::vector<Tri> tris(static_cast<std::size_t>(tri_off[static_cast<std::size_t>(L + 1)]));
    std::vector<Rect> rects(static_cast<std::size_t>(rect_off[static_cast<std::size_t>(L + 1)]));

    std::vector<u32> tri_cursor = tri_off;
    std::vector<u32> rect_cursor = rect_off;
    for (const auto& t : tris_raw) {
        if (t.a < 1 || t.a > L) continue;
        u32& pos = tri_cursor[static_cast<std::size_t>(t.a)];
        tris[static_cast<std::size_t>(pos++)] = t;
    }
    for (const auto& r : rects_raw) {
        if (r.a < 1 || r.a > L) continue;
        u32& pos = rect_cursor[static_cast<std::size_t>(r.a)];
        rects[static_cast<std::size_t>(pos++)] = r;
    }

    std::vector<GroupRange> groups;
    groups.reserve(std::min(tris.size(), rects.size()) / 16);
    for (int a = 1; a <= L; ++a) {
        const std::size_t tb = static_cast<std::size_t>(tri_off[static_cast<std::size_t>(a)]);
        const std::size_t te = static_cast<std::size_t>(tri_off[static_cast<std::size_t>(a + 1)]);
        const std::size_t rb = static_cast<std::size_t>(rect_off[static_cast<std::size_t>(a)]);
        const std::size_t re = static_cast<std::size_t>(rect_off[static_cast<std::size_t>(a + 1)]);
        if (tb == te || rb == re) continue;
        groups.push_back(GroupRange{a, tb, te, rb, re});
    }
    if (groups.empty()) return 0ULL;

    for (const auto& g : groups) {
        std::sort(tris.begin() + static_cast<std::ptrdiff_t>(g.tri_begin),
                  tris.begin() + static_cast<std::ptrdiff_t>(g.tri_end),
                  [](const Tri& x, const Tri& y) { return x.t < y.t; });
        std::sort(rects.begin() + static_cast<std::ptrdiff_t>(g.rect_begin),
                  rects.begin() + static_cast<std::ptrdiff_t>(g.rect_end),
                  [](const Rect& x, const Rect& y) { return x.h < y.h; });
    }

    if (L <= 100'000) {
        u64 sum = 0ULL;
        std::vector<u32> mark(static_cast<std::size_t>(L + 1), 0U);
        u32 stamp = 0U;
        for (const auto& g : groups) {
            sum += compute_group_sum(g, tris, rects, L, mark, stamp);
        }
        return sum;
    }

    const int threads = detect_thread_count(groups.size());
    if (threads <= 1) {
        u64 sum = 0ULL;
        std::vector<u32> mark(static_cast<std::size_t>(L + 1), 0U);
        u32 stamp = 0U;
        for (const auto& g : groups) {
            sum += compute_group_sum(g, tris, rects, L, mark, stamp);
        }
        return sum;
    }

    std::vector<pthread_t> handles(static_cast<std::size_t>(threads));
    std::vector<WorkerTask> tasks(static_cast<std::size_t>(threads));
    std::atomic<std::size_t> next_idx{0};

    for (int t = 0; t < threads; ++t) {
        auto& task = tasks[static_cast<std::size_t>(t)];
        task.tris = &tris;
        task.rects = &rects;
        task.groups = &groups;
        task.L = L;
        task.next_idx = &next_idx;
        task.mark.clear();
        task.stamp = 0U;
        task.partial = 0ULL;
        pthread_create(&handles[static_cast<std::size_t>(t)], nullptr, worker_entry, &task);
    }

    u64 sum = 0ULL;
    for (int t = 0; t < threads; ++t) {
        pthread_join(handles[static_cast<std::size_t>(t)], nullptr);
        sum += tasks[static_cast<std::size_t>(t)].partial;
    }
    return sum;
}

int main() {
    const u64 s1e4 = S(10000);
    if (s1e4 != 884680ULL) {
        std::cerr << "Validation failed: S(1e4) got " << s1e4 << "\n";
        return 1;
    }

    std::cout << S(10000000ULL) << "\n";
    return 0;
}
