#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using namespace std;

static constexpr uint64_t MOD = 1031ULL * 1031ULL * 1031ULL + 2ULL;

static inline __int128 i128_abs(__int128 x) { return x < 0 ? -x : x; }

static uint64_t isqrt_u64(uint64_t n) {
    long double xd = sqrt((long double)n);
    uint64_t x = (uint64_t)xd;
    while ((__int128)(x + 1) * (x + 1) <= n) ++x;
    while ((__int128)x * x > n) --x;
    return x;
}

static uint64_t icbrt_u64(uint64_t n) {
    long double xd = cbrt((long double)n);
    uint64_t x = (uint64_t)xd;
    auto cube = [](__int128 a) { return a * a * a; };
    while (cube((__int128)(x + 1)) <= n) ++x;
    while (cube((__int128)x) > n) --x;
    return x;
}

static uint64_t i4rt_u64(uint64_t n) {
    long double xd = pow((long double)n, 0.25L);
    uint64_t x = (uint64_t)xd;
    auto pow4 = [](__int128 a) { return a * a * a * a; };
    while (pow4((__int128)(x + 1)) <= n) ++x;
    while (pow4((__int128)x) > n) --x;
    return x;
}

static vector<uint32_t> precompute_cube_free(uint32_t limit) {
    vector<uint32_t> spf(limit + 1);
    for (uint32_t i = 0; i <= limit; ++i) spf[i] = i;
    for (uint32_t i = 2; (uint64_t)i * i <= limit; ++i) {
        if (spf[i] == i) {
            for (uint64_t j = (uint64_t)i * i; j <= limit; j += i) {
                if (spf[(uint32_t)j] == (uint32_t)j) spf[(uint32_t)j] = i;
            }
        }
    }

    vector<uint32_t> cf(limit + 1);
    cf[0] = 0;
    if (limit >= 1) cf[1] = 1;
    for (uint32_t n = 2; n <= limit; ++n) {
        uint32_t x = n;
        uint64_t res = 1;
        while (x > 1) {
            uint32_t p = spf[x];
            uint32_t e = 0;
            while (x % p == 0) {
                x /= p;
                ++e;
            }
            uint32_t r = e % 3;
            if (r == 1) res *= p;
            else if (r == 2) res *= (uint64_t)p * p;
        }
        cf[n] = (uint32_t)res;
    }
    return cf;
}

static inline uint64_t sumsq_mod(uint64_t n) {
    // n(n+1)(2n+1)/6
    __int128 t = (__int128)n * (n + 1) * (2 * (__int128)n + 1);
    t /= 6;
    t %= MOD;
    return (uint64_t)t;
}

struct PairLL {
    long long x;
    long long y;
};

static PairLL example_pair_from_pqm(uint64_t p, uint64_t q, uint64_t m) {
    __int128 X0 = 0;
    __int128 Y0 = 0;
    if (q & 1ULL) {
        __int128 t = (__int128)q + 4 * (__int128)p;
        X0 = (__int128)q * t * t * t;
        __int128 u = (__int128)p - 2 * (__int128)q;
        Y0 = 4 * (__int128)p * u * u * u;
    } else {
        __int128 halfq = (__int128)q / 2;
        __int128 t = halfq + 2 * (__int128)p;
        X0 = 2 * (__int128)q * t * t * t;
        __int128 u = (__int128)p - 2 * (__int128)q;
        Y0 = (__int128)p * u * u * u;
    }
    __int128 mm = (__int128)m * m;
    __int128 X = X0 * mm;
    __int128 Y = Y0 * mm;
    if (i128_abs(X) > i128_abs(Y)) swap(X, Y);
    return {(long long)X, (long long)Y};
}

static uint64_t solve_H_mod(uint64_t N, int threads) {
    uint64_t q_max = i4rt_u64((uint64_t)(4 * (__int128)N));

    uint64_t max_p_even2 = 0;
    if (N >= 4) {
        uint64_t c = icbrt_u64(N / 4);
        if (c > 1) max_p_even2 = (c - 1) / 2;
    }
    uint64_t max_p_odd1 = 0;
    {
        uint64_t c = icbrt_u64(N);
        if (c > 1) max_p_odd1 = (c - 1) / 4;
    }
    uint64_t max_p = max(max_p_even2, max_p_odd1);

    uint32_t limit = (uint32_t)(max<uint64_t>(4 * max_p, 2 * q_max) + 16);
    vector<uint32_t> cubeFree = precompute_cube_free(limit);

    threads = max(1, threads);
    threads = min(threads, 64);

    atomic<uint32_t> nextQ(1);
    vector<uint64_t> partial(threads, 0);

    auto worker = [&](int tid) {
        uint64_t local = 0;
        while (true) {
            uint32_t q = nextQ.fetch_add(1, memory_order_relaxed);
            if (q > q_max) break;

            if (q & 1U) {
                uint64_t Ndiv = N / q;
                uint64_t c = icbrt_u64(Ndiv);
                if (c <= q) continue;
                uint64_t p_max = (c - q) / 4;
                if (p_max == 0) continue;

                uint32_t dx = cubeFree[q];

                for (uint64_t p = 1; p <= p_max; ++p) {
                    if (std::gcd(p, (uint64_t)q) != 1) continue;
                    if (p == 2 * (uint64_t)q) continue;

                    uint64_t t = (uint64_t)q + 4ULL * p;
                    __int128 t3 = (__int128)t * t * t;
                    __int128 X0 = (__int128)q * t3;

                    __int128 u = (__int128)p - 2 * (__int128)q;
                    __int128 Y0 = 4 * (__int128)p * u * u * u;

                    __int128 absY = i128_abs(Y0);
                    __int128 max0 = (X0 > absY) ? X0 : absY;
                    if (max0 > (__int128)N) continue;

                    uint32_t dy = cubeFree[4 * p];
                    if (dx == dy) continue;

                    uint64_t max0_u = (uint64_t)max0;
                    uint64_t mmax = isqrt_u64(N / max0_u);
                    uint64_t s2 = sumsq_mod(mmax);

                    uint64_t absSum = ((uint64_t)X0 + (uint64_t)absY) % MOD;
                    local = (local + (uint64_t)((__int128)absSum * s2 % MOD)) % MOD;
                }
            } else {
                uint64_t Ndiv = N / (2ULL * q);
                uint64_t c = icbrt_u64(Ndiv);
                uint64_t halfq = q / 2ULL;
                if (c <= halfq) continue;
                uint64_t p_max = (c - halfq) / 2;
                if (p_max == 0) continue;

                uint32_t dx = cubeFree[2 * q];

                for (uint64_t p = 1; p <= p_max; ++p) {
                    if (std::gcd(p, (uint64_t)q) != 1) continue;
                    if (p == 2 * (uint64_t)q) continue;

                    uint64_t t = halfq + 2ULL * p;
                    __int128 t3 = (__int128)t * t * t;
                    __int128 X0 = 2 * (__int128)q * t3;

                    __int128 u = (__int128)p - 2 * (__int128)q;
                    __int128 Y0 = (__int128)p * u * u * u;

                    __int128 absY = i128_abs(Y0);
                    __int128 max0 = (X0 > absY) ? X0 : absY;
                    if (max0 > (__int128)N) continue;

                    uint32_t dy = cubeFree[p];
                    if (dx == dy) continue;

                    uint64_t max0_u = (uint64_t)max0;
                    uint64_t mmax = isqrt_u64(N / max0_u);
                    uint64_t s2 = sumsq_mod(mmax);

                    uint64_t absSum = ((uint64_t)X0 + (uint64_t)absY) % MOD;
                    local = (local + (uint64_t)((__int128)absSum * s2 % MOD)) % MOD;
                }
            }
        }
        partial[tid] = local;
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    uint64_t ans = 0;
    for (uint64_t v : partial) ans = (ans + v) % MOD;
    return ans;
}

static void run_selftest() {
    {
        PairLL p1 = example_pair_from_pqm(1, 1, 1);
        assert(p1.x == -4 && p1.y == 125);
    }
    {
        PairLL p2 = example_pair_from_pqm(5, 2, 1);
        assert(p2.x == 5 && p2.y == 5324);
    }
    {
        uint64_t got = solve_H_mod(1000, 1);
        assert(got == 2535);
    }
}

int main(int argc, char** argv) {
    uint64_t N = 1000000000000000ULL;
    int threads = (int)std::thread::hardware_concurrency();
    if (threads <= 0) threads = 4;
    bool selftest = false;

    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        if (s == "--selftest") {
            selftest = true;
        } else if (s.rfind("--threads=", 0) == 0) {
            threads = stoi(s.substr(10));
            threads = max(1, threads);
        } else {
            N = stoull(s);
        }
    }

    if (selftest) run_selftest();

    uint64_t ans = solve_H_mod(N, threads);
    cout << ans << "\n";
    return 0;
}
