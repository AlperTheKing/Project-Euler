#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

std::vector<int> fibonacci_up_to(int limit) {
    std::vector<int> fib{1, 2};
    while (fib.back() < limit) {
        fib.push_back(fib[fib.size() - 1] + fib[fib.size() - 2]);
    }
    return fib;
}

std::vector<int> brute_path(int n) {
    const std::vector<int> fib = fibonacci_up_to(2 * n + 5);
    std::vector<char> is_fib(static_cast<std::size_t>(2 * n + 6), 0);
    for (int x : fib) {
        if (x < static_cast<int>(is_fib.size())) {
            is_fib[static_cast<std::size_t>(x)] = 1;
        }
    }

    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n + 1));
    for (int i = 1; i <= n; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            if (is_fib[static_cast<std::size_t>(i + j)] != 0) {
                adj[static_cast<std::size_t>(i)].push_back(j);
                adj[static_cast<std::size_t>(j)].push_back(i);
            }
        }
    }

    std::vector<int> degree(static_cast<std::size_t>(n + 1), 0);
    for (int i = 1; i <= n; ++i) {
        degree[static_cast<std::size_t>(i)] =
            static_cast<int>(adj[static_cast<std::size_t>(i)].size());
    }

    for (int i = 1; i <= n; ++i) {
        auto& v = adj[static_cast<std::size_t>(i)];
        std::sort(v.begin(), v.end(), [&](int a, int b) {
            if (degree[static_cast<std::size_t>(a)] != degree[static_cast<std::size_t>(b)]) {
                return degree[static_cast<std::size_t>(a)] < degree[static_cast<std::size_t>(b)];
            }
            return a < b;
        });
    }

    std::vector<int> starts;
    for (int i = 1; i <= n; ++i) starts.push_back(i);
    std::sort(starts.begin(), starts.end(), [&](int a, int b) {
        if (degree[static_cast<std::size_t>(a)] != degree[static_cast<std::size_t>(b)]) {
            return degree[static_cast<std::size_t>(a)] < degree[static_cast<std::size_t>(b)];
        }
        return a < b;
    });

    std::vector<char> used(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> path;
    path.reserve(static_cast<std::size_t>(n));

    auto dfs = [&](auto&& self, int v) -> bool {
        if (static_cast<int>(path.size()) == n) {
            return true;
        }
        std::vector<int> cand;
        for (int u : adj[static_cast<std::size_t>(v)]) {
            if (used[static_cast<std::size_t>(u)] == 0) {
                cand.push_back(u);
            }
        }
        std::sort(cand.begin(), cand.end(), [&](int a, int b) {
            int ra = 0;
            int rb = 0;
            for (int x : adj[static_cast<std::size_t>(a)]) {
                if (used[static_cast<std::size_t>(x)] == 0) ++ra;
            }
            for (int x : adj[static_cast<std::size_t>(b)]) {
                if (used[static_cast<std::size_t>(x)] == 0) ++rb;
            }
            if (ra != rb) return ra < rb;
            return a < b;
        });

        for (int u : cand) {
            used[static_cast<std::size_t>(u)] = 1;
            path.push_back(u);
            if (self(self, u)) {
                return true;
            }
            path.pop_back();
            used[static_cast<std::size_t>(u)] = 0;
        }
        return false;
    };

    for (int st : starts) {
        std::fill(used.begin(), used.end(), 0);
        path.clear();
        used[static_cast<std::size_t>(st)] = 1;
        path.push_back(st);
        if (dfs(dfs, st)) {
            if (path.front() > path.back()) {
                std::reverse(path.begin(), path.end());
            }
            return path;
        }
    }

    return {};
}

u64 knight_from_right(u64 right_pos, u64 fk, u64 fk_minus_1) {
    if (right_pos == 1ULL) {
        return fk;
    }
    const u64 m = right_pos / 2ULL;
    const u64 t = static_cast<u64>((static_cast<u128>(m) * static_cast<u128>(fk_minus_1)) % fk);
    u64 v = 0;
    if ((right_pos & 1ULL) == 0ULL) {
        v = t;
    } else {
        v = (fk - t) % fk;
    }
    if (v == 0ULL) {
        v = fk;
    }
    return v;
}

}  // namespace

int main() {
    {
        const std::vector<int> p7 = brute_path(7);
        assert(!p7.empty());
        assert(p7[2] == 7);
    }

    {
        const std::vector<int> p34 = brute_path(34);
        assert(!p34.empty());
        assert(p34[2] == 30);
        for (u64 left = 1; left <= 34; ++left) {
            const u64 right = 34ULL - left + 1ULL;
            const u64 got = knight_from_right(right, 34ULL, 21ULL);
            assert(got == static_cast<u64>(p34[static_cast<std::size_t>(left - 1ULL)]));
        }
    }

    std::vector<u64> fib(85, 0ULL);
    fib[1] = 1ULL;
    fib[2] = 1ULL;
    for (int i = 3; i <= 84; ++i) {
        fib[static_cast<std::size_t>(i)] = fib[static_cast<std::size_t>(i - 1)] + fib[static_cast<std::size_t>(i - 2)];
    }

    const u64 n = fib[83];
    const u64 left_pos = 10'000'000'000'000'000ULL;
    const u64 right_pos = n - left_pos + 1ULL;

    const u64 ans = knight_from_right(right_pos, n, fib[82]);
    std::cout << ans << "\n";
    return 0;
}
