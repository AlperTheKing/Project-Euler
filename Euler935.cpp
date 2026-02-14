#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

static std::vector<u32> build_spf(u32 n) {
    std::vector<u32> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    if (n >= 1) spf[1] = 1;
    for (u32 i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (u32 p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > n) break;
            spf[static_cast<std::size_t>(v)] = p;
            if (p == spf[static_cast<std::size_t>(i)]) break;
        }
    }
    return spf;
}

static int factor_distinct(u32 x, const std::vector<u32>& spf, u32 out[10]) {
    int cnt = 0;
    while (x > 1) {
        const u32 p = spf[static_cast<std::size_t>(x)];
        out[cnt++] = p;
        do {
            x /= p;
        } while (x > 1 && x % p == 0);
    }
    return cnt;
}

static u64 count_relprimes_pie(u32 pmax, u32 k, const std::vector<u32>& spf) {
    u32 fac[10];
    const int m = factor_distinct(k, spf, fac);
    const int lim = 1 << m;
    u64 res = 0;
    for (int mask = 0; mask < lim; ++mask) {
        u64 prod = 1;
        int bits = 0;
        int mm = mask;
        bool ok = true;
        while (mm) {
            const int b = __builtin_ctz(mm);
            mm &= (mm - 1);
            ++bits;
            const u32 p = fac[b];
            if (prod > pmax / p) {
                ok = false;
                break;
            }
            prod *= p;
        }
        if (!ok) continue;
        const u64 term = pmax / prod;
        if (bits & 1) {
            res -= term;
        } else {
            res += term;
        }
    }
    return res;
}

static u64 solve(u32 n_max) {
    const u32 k_max = n_max + 4;
    const std::vector<u32> spf = build_spf(k_max);
    const u32 q_max = (k_max - 1) / 4;

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;
    if (threads > 12) threads = 12;
    if (q_max < 200'000) threads = 1;

    std::vector<std::thread> workers;
    std::vector<u64> partial(static_cast<std::size_t>(threads), 0);
    workers.reserve(threads);

    const u32 chunk = (q_max + threads - 1) / threads;
    for (unsigned t = 0; t < threads; ++t) {
        const u32 lo = static_cast<u32>(t) * chunk + 1;
        const u32 hi = std::min<u32>(q_max, static_cast<u32>(t + 1) * chunk);
        if (lo > hi) continue;
        workers.emplace_back([&, t, lo, hi]() {
            u64 local = 0;
            for (u32 q = lo; q <= hi; ++q) {
                u32 k = (q << 2) + 1;
                u32 n = k - 1;
                if (n <= n_max) local += count_relprimes_pie(q, k, spf);

                k = (q << 2) + 2;
                n = k - 2;
                if (n <= n_max) local += count_relprimes_pie(q, k, spf);

                k = (q << 2) + 3;
                n = k - 1;
                if (n <= n_max) local += count_relprimes_pie(q, k, spf);

                k = (q << 2) + 4;
                n = k - 4;
                if (n <= n_max) local += count_relprimes_pie(q, k, spf);
            }
            partial[static_cast<std::size_t>(t)] = local;
        });
    }

    for (auto& th : workers) th.join();
    u64 ans = 0;
    for (u64 v : partial) ans += v;
    return ans;
}

int main() {
    assert(solve(6) == 4);
    assert(solve(100) == 805);
    std::cout << solve(100'000'000) << '\n';
    return 0;
}
