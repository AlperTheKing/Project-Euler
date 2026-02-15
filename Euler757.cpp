#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

u64 isqrt_u64(const u64 n) {
    long double d = static_cast<long double>(n);
    u64 x = static_cast<u64>(std::sqrt(d));
    while ((x + 1ULL) * (x + 1ULL) <= n) {
        ++x;
    }
    while (x * x > n) {
        --x;
    }
    return x;
}

u64 max_y_with_pronic_leq(const u64 t) {
    const u64 s = isqrt_u64(1ULL + 4ULL * t);
    return (s - 1ULL) / 2ULL;
}

u64 count_stealthy(const u64 limit) {
    const u64 y_global_max = max_y_with_pronic_leq(limit);
    std::vector<u64> pronic(static_cast<std::size_t>(y_global_max + 1ULL), 0ULL);
    for (u64 y = 1ULL; y <= y_global_max; ++y) {
        pronic[static_cast<std::size_t>(y)] = y * (y + 1ULL);
    }

    u64 x_max = 0ULL;
    for (u64 x = 1ULL; x <= y_global_max; ++x) {
        const u64 p = pronic[static_cast<std::size_t>(x)];
        if (p > limit / p) {
            break;
        }
        x_max = x;
    }

    struct Node {
        u64 value;
        u64 px;
        u32 y;
        u32 y_max;
    };

    std::vector<Node> heap;
    heap.reserve(static_cast<std::size_t>(x_max));

    auto sift_up = [&](std::size_t idx) {
        while (idx > 0U) {
            const std::size_t parent = (idx - 1U) >> 1U;
            if (heap[parent].value <= heap[idx].value) {
                break;
            }
            std::swap(heap[parent], heap[idx]);
            idx = parent;
        }
    };

    auto sift_down = [&](std::size_t idx) {
        const std::size_t n = heap.size();
        while (true) {
            std::size_t left = (idx << 1U) + 1U;
            if (left >= n) {
                break;
            }
            std::size_t right = left + 1U;
            std::size_t best = left;
            if (right < n && heap[right].value < heap[left].value) {
                best = right;
            }
            if (heap[idx].value <= heap[best].value) {
                break;
            }
            std::swap(heap[idx], heap[best]);
            idx = best;
        }
    };

    for (u64 x = 1ULL; x <= x_max; ++x) {
        const u64 px = pronic[static_cast<std::size_t>(x)];
        const u64 max_t = limit / px;
        const u64 y_lim = max_y_with_pronic_leq(max_t);
        if (y_lim < x) {
            continue;
        }
        heap.push_back(Node{px * px, px, static_cast<u32>(x), static_cast<u32>(y_lim)});
        sift_up(heap.size() - 1U);
    }

    u64 count = 0ULL;
    u64 last = static_cast<u64>(-1);

    while (!heap.empty()) {
        Node cur = heap[0];

        if (cur.value != last) {
            ++count;
            last = cur.value;
        }

        if (cur.y < cur.y_max) {
            ++cur.y;
            cur.value = cur.px * pronic[static_cast<std::size_t>(cur.y)];
            heap[0] = cur;
            sift_down(0U);
        } else {
            heap[0] = heap.back();
            heap.pop_back();
            if (!heap.empty()) {
                sift_down(0U);
            }
        }
    }

    return count;
}

}  // namespace

int main() {
    assert(count_stealthy(1'000'000ULL) == 2'851ULL);
    std::cout << count_stealthy(100'000'000'000'000ULL) << '\n';
    return 0;
}
