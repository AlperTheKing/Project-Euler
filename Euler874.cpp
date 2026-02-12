#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

using i64 = long long;
using u64 = std::uint64_t;

static std::vector<int> first_primes(int count) {
    if (count <= 0) return {};
    int limit;
    if (count < 6) {
        limit = 15;
    } else {
        double n = static_cast<double>(count);
        limit = static_cast<int>(n * (std::log(n) + std::log(std::log(n)))) + 32;
    }

    while (true) {
        std::vector<bool> is_prime(limit + 1, true);
        is_prime[0] = false;
        is_prime[1] = false;
        for (int i = 2; 1LL * i * i <= limit; ++i) {
            if (!is_prime[i]) continue;
            for (int j = i * i; j <= limit; j += i) is_prime[j] = false;
        }

        std::vector<int> primes;
        primes.reserve(count);
        for (int x = 2; x <= limit && static_cast<int>(primes.size()) < count; ++x) {
            if (is_prime[x]) primes.push_back(x);
        }
        if (static_cast<int>(primes.size()) >= count) return primes;
        limit *= 2;
    }
}

static u64 maximal_score_fast(int k, u64 n, const std::vector<int>& primes) {
    assert(k >= 2);
    assert(static_cast<int>(primes.size()) >= k);
    assert(n >= static_cast<u64>(k - 1));

    int top = primes[k - 1];
    std::vector<int> loss(k, 0);
    for (int d = 1; d < k; ++d) {
        loss[d] = top - primes[k - 1 - d];
    }

    int residue = static_cast<int>((k - (n % static_cast<u64>(k))) % static_cast<u64>(k));

    const u64 INF = std::numeric_limits<u64>::max() / 4;
    std::vector<u64> dist(k, INF);
    std::priority_queue<std::pair<u64, int>, std::vector<std::pair<u64, int>>, std::greater<>> pq;

    dist[0] = 0;
    pq.push({0, 0});

    while (!pq.empty()) {
        auto [du, u] = pq.top();
        pq.pop();
        if (du != dist[u]) continue;

        for (int d = 1; d < k; ++d) {
            int v = u + d;
            if (v >= k) v -= k;
            u64 nd = du + static_cast<u64>(loss[d]);
            if (nd < dist[v]) {
                dist[v] = nd;
                pq.push({nd, v});
            }
        }
    }

    return n * static_cast<u64>(top) - dist[residue];
}

static u64 maximal_score_dp(int k, int n, const std::vector<int>& primes) {
    const i64 NEG = std::numeric_limits<i64>::min() / 4;
    std::vector<i64> dp(k, NEG), ndp(k, NEG);
    dp[0] = 0;

    for (int i = 0; i < n; ++i) {
        std::fill(ndp.begin(), ndp.end(), NEG);
        for (int r = 0; r < k; ++r) {
            if (dp[r] == NEG) continue;
            for (int a = 0; a < k; ++a) {
                int nr = (r + a) % k;
                ndp[nr] = std::max(ndp[nr], dp[r] + primes[a]);
            }
        }
        dp.swap(ndp);
    }

    return static_cast<u64>(dp[0]);
}

int main() {
    {
        auto p = first_primes(10);
        assert(maximal_score_fast(2, 5, p) == 14ULL);
        for (int k = 2; k <= 8; ++k) {
            for (int n = k - 1; n <= k + 5; ++n) {
                assert(maximal_score_fast(k, static_cast<u64>(n), p) == maximal_score_dp(k, n, p));
            }
        }
    }

    constexpr int k = 7000;
    auto primes = first_primes(k + 1);
    u64 n = static_cast<u64>(primes[k]);
    std::cout << maximal_score_fast(k, n, primes) << '\n';
    return 0;
}
