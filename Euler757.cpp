#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <queue>
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

    std::vector<u64> y_max(static_cast<std::size_t>(x_max + 1ULL), 0ULL);

    struct Node {
        u64 value;
        u32 x;
        u32 y;
    };
    auto cmp = [](const Node& lhs, const Node& rhs) { return lhs.value > rhs.value; };
    std::priority_queue<Node, std::vector<Node>, decltype(cmp)> pq(cmp);

    for (u64 x = 1ULL; x <= x_max; ++x) {
        const u64 max_t = limit / pronic[static_cast<std::size_t>(x)];
        const u64 y_lim = max_y_with_pronic_leq(max_t);
        if (y_lim < x) {
            continue;
        }
        y_max[static_cast<std::size_t>(x)] = y_lim;
        const u64 value = pronic[static_cast<std::size_t>(x)] * pronic[static_cast<std::size_t>(x)];
        pq.push(Node{value, static_cast<u32>(x), static_cast<u32>(x)});
    }

    u64 count = 0ULL;
    u64 last = 0ULL;
    bool has_last = false;

    while (!pq.empty()) {
        const Node cur = pq.top();
        pq.pop();

        if (!has_last || cur.value != last) {
            ++count;
            last = cur.value;
            has_last = true;
        }

        const u64 x = static_cast<u64>(cur.x);
        const u64 y_next = static_cast<u64>(cur.y) + 1ULL;
        if (y_next <= y_max[static_cast<std::size_t>(x)]) {
            const u64 value =
                pronic[static_cast<std::size_t>(x)] * pronic[static_cast<std::size_t>(y_next)];
            pq.push(Node{value, cur.x, static_cast<u32>(y_next)});
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
