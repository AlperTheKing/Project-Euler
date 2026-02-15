#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;

namespace {

constexpr u64 MOD = 1000000007ULL;

inline void add_mod(u64& x, u64 y) {
    x += y;
    if (x >= MOD) x -= MOD;
}

u64 compute_f(u64 n, int m, bool require_stable) {
    const int dim_inv = m + 1;
    const int dim_last132 = m + 3;
    const int dim_last21 = m + 3;
    const std::size_t plane = static_cast<std::size_t>(dim_last132) * static_cast<std::size_t>(dim_last21);
    const std::size_t total_states = static_cast<std::size_t>(dim_inv) * plane;

    auto idx = [=](int inv, int last132, int last21) -> std::size_t {
        return static_cast<std::size_t>(inv) * plane +
               static_cast<std::size_t>(last132 + 1) * static_cast<std::size_t>(dim_last21) +
               static_cast<std::size_t>(last21);
    };

    std::vector<u64> cur(total_states, 0);
    std::vector<u64> nxt(total_states, 0);
    cur[idx(m, -1, 0)] = 1;

    u64 total = 1;
    const int max_len = m + 2;
    const int limit = (n < static_cast<u64>(max_len)) ? static_cast<int>(n) : max_len;

    for (int length = 1; length <= limit; ++length) {
        std::fill(nxt.begin(), nxt.end(), 0);
        for (int inv = 0; inv <= m; ++inv) {
            const int upper = (inv + 1 < length) ? (inv + 1) : length;
            if (upper <= 0) continue;
            for (int last132p1 = 0; last132p1 <= length; ++last132p1) {
                const int start = last132p1;
                if (start >= upper) continue;
                for (int last21 = 0; last21 < length; ++last21) {
                    const u64 cnt = cur[static_cast<std::size_t>(inv) * plane +
                                        static_cast<std::size_t>(last132p1) * static_cast<std::size_t>(dim_last21) +
                                        static_cast<std::size_t>(last21)];
                    if (cnt == 0) continue;
                    for (int i = start; i < upper; ++i) {
                        add_mod(total, cnt % MOD);
                        if (i < last21) {
                            add_mod(nxt[idx(inv - i, i, last21 + 1)], cnt % MOD);
                        } else {
                            add_mod(nxt[idx(inv - i, last132p1 - 1, i)], cnt % MOD);
                        }
                    }
                }
            }
        }
        cur.swap(nxt);
    }

    if (n <= static_cast<u64>(max_len)) return total;

    u64 g = 0;
    for (u64 v : cur) add_mod(g, v % MOD);

    if (require_stable) {
        const int length = max_len + 1;
        std::fill(nxt.begin(), nxt.end(), 0);
        for (int inv = 0; inv <= m; ++inv) {
            const int upper = (inv + 1 < length) ? (inv + 1) : length;
            if (upper <= 0) continue;
            for (int last132p1 = 0; last132p1 <= max_len; ++last132p1) {
                const int start = last132p1;
                if (start >= upper) continue;
                for (int last21 = 0; last21 <= max_len; ++last21) {
                    const u64 cnt = cur[static_cast<std::size_t>(inv) * plane +
                                        static_cast<std::size_t>(last132p1) * static_cast<std::size_t>(dim_last21) +
                                        static_cast<std::size_t>(last21)];
                    if (cnt == 0) continue;
                    for (int i = start; i < upper; ++i) {
                        if (i < last21) {
                            add_mod(nxt[idx(inv - i, i, last21 + 1)], cnt % MOD);
                        } else {
                            add_mod(nxt[idx(inv - i, last132p1 - 1, i)], cnt % MOD);
                        }
                    }
                }
            }
        }
        u64 g_next = 0;
        for (u64 v : nxt) add_mod(g_next, v % MOD);
        if (g_next != g) {
            std::cerr << "Stabilization check failed.\n";
            std::exit(1);
        }
    }

    const u64 extra = ((n - static_cast<u64>(max_len)) % MOD) * g % MOD;
    return (total + extra) % MOD;
}

void run_checks() {
    struct Check {
        u64 n;
        int m;
        u64 expected;
    };
    const Check checks[] = {
        {2, 0, 3},
        {4, 5, 32},
        {10, 25, 294400},
    };
    for (const auto& chk : checks) {
        const u64 got = compute_f(chk.n, chk.m, false);
        if (got != chk.expected) {
            std::cerr << "Validation failed for f(" << chk.n << "," << chk.m << "): "
                      << got << " != " << chk.expected << "\n";
            std::exit(1);
        }
    }
}

}  // namespace

int main() {
    run_checks();
    const u64 n = 1000000000000000000ULL;
    const int m = 40;
    std::cout << compute_f(n, m, true) << '\n';
    return 0;
}
