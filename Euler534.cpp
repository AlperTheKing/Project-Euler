#include <cstdint>
#include <iomanip>
#include <iostream>
#include <atomic>
#include <algorithm>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

u64 count_configurations(int n, int D) {
    // Place exactly one queen per row. A queen attacks vertically/diagonally only up to D rows away.
    // Row-by-row DP with state = last D columns (4 bits each, most recent at LSB).
    const u64 all_cols = (n == 64) ? ~0ULL : ((1ULL << n) - 1ULL);
    const u64 mask = (D == 0) ? 0ULL : ((1ULL << (4 * D)) - 1ULL);

    struct Node {
        u64 state;
        u64 ways;
    };

    std::vector<Node> dp, next, generated;
    dp.reserve(1 << 18);
    next.reserve(1 << 18);
    generated.reserve(1 << 20);
    dp.push_back(Node{0ULL, 1ULL});

    for (int row = 0; row < n; ++row) {
        generated.clear();
        const int len = (row < D) ? row : D;
        for (const Node& node : dp) {
            const u64 state = node.state;
            const u64 ways = node.ways;
            u64 blocked = 0;
            for (int dist = 1; dist <= len; ++dist) {
                const int c_prev = static_cast<int>((state >> (4 * (dist - 1))) & 0xFULL);
                blocked |= 1ULL << c_prev;
                const int c1 = c_prev + dist;
                const int c2 = c_prev - dist;
                if (c1 >= 0 && c1 < n) {
                    blocked |= 1ULL << c1;
                }
                if (c2 >= 0 && c2 < n) {
                    blocked |= 1ULL << c2;
                }
            }
            u64 avail = all_cols & ~blocked;
            while (avail) {
                const int col = __builtin_ctzll(avail);
                avail &= (avail - 1);
                const u64 ns = (((state << 4) | static_cast<u64>(col)) & mask);
                generated.push_back(Node{ns, ways});
            }
        }

        std::sort(generated.begin(), generated.end(),
                  [](const Node& a, const Node& b) { return a.state < b.state; });
        next.clear();
        next.reserve(generated.size());
        for (const Node& g : generated) {
            if (!next.empty() && next.back().state == g.state) {
                next.back().ways += g.ways;
            } else {
                next.push_back(g);
            }
        }

        dp.swap(next);
    }

    u128 total = 0;
    for (const Node& node : dp) {
        total += static_cast<u128>(node.ways);
    }
    return static_cast<u64>(total);
}

u64 S(int n) {
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 8;
    }
    threads = std::min<unsigned>(threads, static_cast<unsigned>(n));

    std::atomic<int> next_w{0};
    std::vector<u128> partial(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            u128 local = 0;
            while (true) {
                const int w = next_w.fetch_add(1, std::memory_order_relaxed);
                if (w >= n) {
                    break;
                }
                const int D = (n - 1) - w;
                local += static_cast<u128>(count_configurations(n, D));
            }
            partial[tid] = local;
        });
    }
    for (auto& t : workers) {
        t.join();
    }

    u128 sum = 0;
    for (u128 v : partial) {
        sum += v;
    }

    return static_cast<u64>(sum);
}

bool run_checkpoints() {
    if (count_configurations(4, 3) != 2ULL) {  // w=0 => D=3
        std::cerr << "Checkpoint failed: Q(4,0)\n";
        return false;
    }
    if (count_configurations(4, 1) != 16ULL) {  // w=2 => D=1
        std::cerr << "Checkpoint failed: Q(4,2)\n";
        return false;
    }
    if (count_configurations(4, 0) != 256ULL) {  // w=3 => D=0
        std::cerr << "Checkpoint failed: Q(4,3)\n";
        return false;
    }
    if (S(4) != 276ULL) {
        std::cerr << "Checkpoint failed: S(4)\n";
        return false;
    }
    if (S(5) != 3347ULL) {
        std::cerr << "Checkpoint failed: S(5)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }
    std::cout << S(14) << '\n';
    return 0;
}
