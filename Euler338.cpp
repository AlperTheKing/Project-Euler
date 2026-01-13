#include <algorithm>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

using int64 = long long;
using i128 = __int128_t;

namespace {

constexpr int64 MOD = 100000000LL;

int64 iroot3(int64 n) {
    long double approx = cbrt(static_cast<long double>(n));
    int64 x = static_cast<int64>(approx);
    auto cube = [](int64 v) -> i128 {
        return static_cast<i128>(v) * v * v;
    };
    while (cube(x + 1) <= n) ++x;
    while (cube(x) > n) --x;
    return x;
}

int64 divisor_summatory(int64 n, std::unordered_map<int64, int64>& cache) {
    if (n <= 0) return 0;
    auto it = cache.find(n);
    if (it != cache.end()) return it->second;

    int64 res = 0;
    int64 k = 1;
    while (k <= n) {
        int64 q = n / k;
        int64 r = n / q;
        res += q * (r - k + 1);
        k = r + 1;
    }

    cache[n] = res;
    return res;
}

i128 sum_floor_products(int64 n) {
    if (n < 2) return 0;
    i128 sum = 0;
    int64 i = 2;
    while (i <= n) {
        int64 a = n / i;
        int64 b = n / (i - 1);

        int64 next_a = n / a + 1;
        int64 next_b = n / b + 2;
        int64 nxt = std::min(next_a, next_b);
        if (nxt > n + 1) nxt = n + 1;

        int64 cnt = nxt - i;
        sum += static_cast<i128>(cnt) * a * b;
        i = nxt;
    }
    return sum;
}

i128 sum_double_floor(int64 n, int64 limit) {
    // sum_{a=1..limit} sum_{b=1..limit} floor(n/(a*b))
    int thread_count = static_cast<int>(std::thread::hardware_concurrency());
    if (thread_count <= 0) thread_count = 1;
    if (limit < thread_count) thread_count = static_cast<int>(limit);

    std::vector<i128> partial(thread_count, 0);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    for (int t = 0; t < thread_count; ++t) {
        int64 start = 1 + (limit * t) / thread_count;
        int64 end = (limit * (t + 1)) / thread_count;
        workers.emplace_back([&, t, start, end]() {
            i128 local = 0;
            for (int64 a = start; a <= end; ++a) {
                for (int64 b = 1; b <= limit; ++b) {
                    int64 denom = a * b;
                    local += n / denom;
                }
            }
            partial[t] = local;
        });
    }
    for (auto& th : workers) th.join();

    i128 total = 0;
    for (auto v : partial) total += v;
    return total;
}

i128 triple_count(int64 n) {
    // T(n) = #{(a,b,c) >= 1 : a*b*c <= n}
    // Use inclusion-exclusion with L = floor(n^(1/3)).
    int64 L = iroot3(n);

    std::unordered_map<int64, int64> cache;
    cache.reserve(static_cast<size_t>(L * 2 + 100));

    i128 sum_a = 0;
    int64 a = 1;
    while (a <= L) {
        int64 q = n / a;
        int64 r = n / q;
        if (r > L) r = L;
        int64 cnt = r - a + 1;
        sum_a += static_cast<i128>(cnt) * divisor_summatory(q, cache);
        a = r + 1;
    }

    i128 sum_ab = sum_double_floor(n, L);
    i128 l3 = static_cast<i128>(L) * L * L;

    return 3 * sum_a - 3 * sum_ab + l3;
}

i128 compute_G(int64 n) {
    if (n < 2) return 0;

    std::unordered_map<int64, int64> cache;
    cache.reserve(1 << 16);

    i128 s_prod = sum_floor_products(n);
    int64 d_n = divisor_summatory(n, cache);
    i128 t_n = triple_count(n);

    return s_prod - t_n + d_n;
}

int64 mod_norm(i128 v) {
    v %= MOD;
    if (v < 0) v += MOD;
    return static_cast<int64>(v);
}

}  // namespace

int main() {
    const int64 N = 1000000000000LL;

    struct Check { int64 n; int64 expected; };
    const Check checks[] = {
        {10, 55},
        {1000, 971745},
        {100000, 9992617687LL},
    };

    for (const auto& chk : checks) {
        int64 got = mod_norm(compute_G(chk.n));
        if (got != (chk.expected % MOD)) {
            std::cerr << "Validation failure: G(" << chk.n << ") mod 1e8 = "
                      << got << ", expected " << (chk.expected % MOD) << '\n';
            return 1;
        }
    }

    int64 answer = mod_norm(compute_G(N));
    std::cout << answer << '\n';
    return 0;
}
