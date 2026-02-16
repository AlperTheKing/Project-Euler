#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;

namespace {

constexpr uint64_t MOD = 100000000ULL;
constexpr long long DEFAULT_N = 100000000000LL;
constexpr long long PRECOMP_LIMIT = 25000000LL;

struct FastHash {
    size_t operator()(uint64_t x) const noexcept {
        static const uint64_t kMul = 0x9e3779b97f4a7c15ULL;
        x += kMul;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<size_t>(x ^ (x >> 31));
    }
};

inline uint64_t mod_norm(int64_t v) {
    v %= static_cast<int64_t>(MOD);
    if (v < 0) v += static_cast<int64_t>(MOD);
    return static_cast<uint64_t>(v);
}

inline uint64_t mod_mul(uint64_t a, uint64_t b) {
    return static_cast<uint64_t>((__int128)(a % MOD) * (b % MOD) % MOD);
}

uint64_t pow2_mod(unsigned __int128 exp) {
    uint64_t base = 2 % MOD;
    uint64_t res = 1;
    while (exp > 0) {
        if (exp & 1) res = mod_mul(res, base);
        base = mod_mul(base, base);
        exp >>= 1;
    }
    return res;
}

uint64_t sum1_upto(long long n) {
    if (n <= 0) return 0;
    __int128 v = static_cast<__int128>(n) * (n + 1) / 2;
    return static_cast<uint64_t>(v % MOD);
}

uint64_t sum2_upto(long long n) {
    if (n <= 0) return 0;
    __int128 v = static_cast<__int128>(n) * (n + 1) * (2 * n + 1) / 6;
    return static_cast<uint64_t>(v % MOD);
}

uint64_t sum3_upto(long long n) {
    if (n <= 0) return 0;
    uint64_t a = static_cast<uint64_t>(n % static_cast<long long>(MOD));
    uint64_t b = static_cast<uint64_t>((n + 1) % static_cast<long long>(MOD));
    if ((n & 1LL) == 0) {
        a = static_cast<uint64_t>((n / 2) % static_cast<long long>(MOD));
    } else {
        b = static_cast<uint64_t>(((n + 1) / 2) % static_cast<long long>(MOD));
    }
    uint64_t t = mod_mul(a, b);
    return mod_mul(t, t);
}

uint64_t sum_range_1(long long l, long long r) {
    return mod_norm(static_cast<int64_t>(sum1_upto(r)) - static_cast<int64_t>(sum1_upto(l - 1)));
}

uint64_t sum_range_2(long long l, long long r) {
    return mod_norm(static_cast<int64_t>(sum2_upto(r)) - static_cast<int64_t>(sum2_upto(l - 1)));
}

uint64_t sum_range_3(long long l, long long r) {
    return mod_norm(static_cast<int64_t>(sum3_upto(r)) - static_cast<int64_t>(sum3_upto(l - 1)));
}

uint64_t sum_q_mod(long long l, long long r) {
    __int128 cnt = static_cast<__int128>(r - l + 1);
    __int128 sum = static_cast<__int128>(l + r) * cnt / 2;
    return static_cast<uint64_t>(sum % MOD);
}

uint64_t sum_lpow_upto(long long n, uint64_t pow2_n_plus1) {
    if (n < 0) return 0;
    uint64_t n_mod = static_cast<uint64_t>(n % static_cast<long long>(MOD));
    uint64_t coef = mod_norm(static_cast<int64_t>(n_mod) - 1);
    return (mod_mul(coef, pow2_n_plus1) + 2) % MOD;
}

uint64_t sum_l2pow_upto(long long n, uint64_t pow2_n_plus1) {
    if (n < 0) return 0;
    uint64_t n_mod = static_cast<uint64_t>(n % static_cast<long long>(MOD));
    uint64_t n2_mod = mod_mul(n_mod, n_mod);
    uint64_t coef = mod_norm(static_cast<int64_t>(n2_mod) -
                             static_cast<int64_t>(mod_mul(2, n_mod)) + 3);
    uint64_t val = mod_mul(coef, pow2_n_plus1);
    return mod_norm(static_cast<int64_t>(val) - 6);
}

struct TotientSummatory {
    int limit = 0;
    vector<int> phi;
    vector<uint32_t> pref_phi;
    vector<uint32_t> pref_weighted;
    vector<uint32_t> pref_weighted2;
    unordered_map<long long, uint64_t, FastHash> cache_phi;
    unordered_map<long long, uint64_t, FastHash> cache_weighted;
    unordered_map<long long, uint64_t, FastHash> cache_weighted2;

    explicit TotientSummatory(long long max_n) {
        limit = static_cast<int>(min(max_n, PRECOMP_LIMIT));
        if (limit < 1) limit = 1;
        phi.assign(limit + 1, 0);
        pref_phi.assign(limit + 1, 0);
        pref_weighted.assign(limit + 1, 0);
        pref_weighted2.assign(limit + 1, 0);

        vector<int> primes;
        vector<uint8_t> is_comp(limit + 1, 0);
        phi[1] = 1;
        for (int i = 2; i <= limit; ++i) {
            if (!is_comp[i]) {
                primes.push_back(i);
                phi[i] = i - 1;
            }
            for (int p : primes) {
                int64_t v = 1LL * i * p;
                if (v > limit) break;
                is_comp[static_cast<int>(v)] = 1;
                if (i % p == 0) {
                    phi[static_cast<int>(v)] = phi[i] * p;
                    break;
                }
                phi[static_cast<int>(v)] = phi[i] * (p - 1);
            }
        }

        for (int i = 1; i <= limit; ++i) {
            uint64_t add_phi = static_cast<uint64_t>(phi[i]) % MOD;
            pref_phi[i] = static_cast<uint32_t>((pref_phi[i - 1] + add_phi) % MOD);

            uint64_t add_weight = (static_cast<uint64_t>(i) * phi[i]) % MOD;
            pref_weighted[i] = static_cast<uint32_t>((pref_weighted[i - 1] + add_weight) % MOD);

            uint64_t i_mod = static_cast<uint64_t>(i) % MOD;
            uint64_t i2_mod = mod_mul(i_mod, i_mod);
            uint64_t add_weight2 = mod_mul(i2_mod, static_cast<uint64_t>(phi[i]) % MOD);
            pref_weighted2[i] = static_cast<uint32_t>((pref_weighted2[i - 1] + add_weight2) % MOD);
        }

        cache_phi.reserve(1 << 20);
        cache_weighted.reserve(1 << 20);
        cache_weighted2.reserve(1 << 20);
    }

    uint64_t phi_sum(long long n) {
        if (n <= limit) return pref_phi[static_cast<int>(n)];
        auto it = cache_phi.find(n);
        if (it != cache_phi.end()) return it->second;

        uint64_t res = sum1_upto(n);
        for (long long l = 2; l <= n; ) {
            long long t = n / l;
            long long r = n / t;
            uint64_t sub = mod_mul(static_cast<uint64_t>(r - l + 1), phi_sum(t));
            res = (res + MOD - sub) % MOD;
            l = r + 1;
        }
        cache_phi[n] = res;
        return res;
    }

    uint64_t phi_weighted_sum(long long n) {
        if (n <= limit) return pref_weighted[static_cast<int>(n)];
        auto it = cache_weighted.find(n);
        if (it != cache_weighted.end()) return it->second;

        uint64_t res = sum2_upto(n);
        for (long long l = 2; l <= n; ) {
            long long t = n / l;
            long long r = n / t;
            uint64_t sum_q = sum_q_mod(l, r);
            uint64_t sub = mod_mul(sum_q, phi_weighted_sum(t));
            res = (res + MOD - sub) % MOD;
            l = r + 1;
        }
        cache_weighted[n] = res;
        return res;
    }

    uint64_t phi_weighted2_sum(long long n) {
        if (n <= limit) return pref_weighted2[static_cast<int>(n)];
        auto it = cache_weighted2.find(n);
        if (it != cache_weighted2.end()) return it->second;

        uint64_t res = sum3_upto(n);
        for (long long l = 2; l <= n; ) {
            long long t = n / l;
            long long r = n / t;
            uint64_t sum_q2 = sum_range_2(l, r);
            uint64_t sub = mod_mul(sum_q2, phi_weighted2_sum(t));
            res = (res + MOD - sub) % MOD;
            l = r + 1;
        }
        cache_weighted2[n] = res;
        return res;
    }
};

struct RawRange {
    long long L1;
    long long L2;
    long long M;
    long long M2;
};

struct Range {
    long long L1;
    long long L2;
    uint64_t c1;
    uint64_t sxy1;
    uint64_t sprod1;
    uint64_t c2;
    uint64_t sxy2;
    uint64_t sprod2;
};

uint64_t compute_T(long long N, TotientSummatory& ts, unsigned threads) {
    unsigned __int128 exp = static_cast<unsigned __int128>(N + 1) * (N + 1);
    uint64_t total_subsets = pow2_mod(exp);

    uint64_t points_mod = mod_mul(static_cast<uint64_t>((N + 1) % MOD),
                                  static_cast<uint64_t>((N + 1) % MOD));

    if (N < 2) {
        uint64_t non_titanic = (1 + points_mod) % MOD;
        return mod_norm(static_cast<int64_t>(total_subsets) - static_cast<int64_t>(non_titanic));
    }

    vector<RawRange> raw_ranges;
    vector<long long> ms;
    for (long long L = 2; L <= N; ) {
        long long M1 = N / L;
        long long M2 = (L + 1 <= N) ? (N / (L + 1)) : 0;
        long long R1 = N / M1;
        long long R2 = (M2 == 0) ? N : (N / M2 - 1);
        long long R = min(R1, R2);
        raw_ranges.push_back({L, R, M1, M2});
        ms.push_back(M1);
        if (M2 > 0) ms.push_back(M2);
        L = R + 1;
    }

    sort(ms.begin(), ms.end());
    ms.erase(unique(ms.begin(), ms.end()), ms.end());

    struct PrecompVals {
        uint64_t c;
        uint64_t sxy;
        uint64_t sprod;
    };

    unordered_map<long long, PrecompVals, FastHash> precomp;
    precomp.reserve(ms.size() * 2);
    for (long long M : ms) {
        uint64_t sum_phi = ts.phi_sum(M);
        uint64_t sum_weighted = ts.phi_weighted_sum(M);
        uint64_t sum_weighted2 = ts.phi_weighted2_sum(M);
        uint64_t c = (2 * sum_phi + MOD - 1) % MOD;
        uint64_t sxy = (3 * sum_weighted + MOD - 1) % MOD;
        uint64_t sprod = sum_weighted2 % MOD;
        precomp.emplace(M, PrecompVals{c, sxy, sprod});
    }

    vector<Range> ranges;
    ranges.reserve(raw_ranges.size());
    for (const auto& rr : raw_ranges) {
        auto it1 = precomp.find(rr.M);
        PrecompVals v1 = it1->second;
        PrecompVals v2{0, 0, 0};
        if (rr.M2 > 0) {
            auto it2 = precomp.find(rr.M2);
            v2 = it2->second;
        }
        ranges.push_back({rr.L1, rr.L2, v1.c, v1.sxy, v1.sprod, v2.c, v2.sxy, v2.sprod});
    }

    uint64_t a = static_cast<uint64_t>((N + 1) % static_cast<long long>(MOD));
    uint64_t a2 = mod_mul(a, a);

    if (threads == 0) threads = 1;
    if (ranges.size() < 2000) threads = 1;
    if (threads > ranges.size()) threads = static_cast<unsigned>(ranges.size());

    vector<uint64_t> partial(threads, 0);
    atomic<size_t> next{0};
    const size_t chunk = 128;

    auto worker = [&](unsigned idx) {
        uint64_t local = 0;
        while (true) {
            size_t start = next.fetch_add(chunk, memory_order_relaxed);
            if (start >= ranges.size()) break;
            size_t end = min(start + chunk, ranges.size());
            for (size_t i = start; i < end; ++i) {
                const Range& rg = ranges[i];
                uint64_t q2 = mod_norm(static_cast<int64_t>(rg.sprod1) - static_cast<int64_t>(rg.sprod2));
                uint64_t q1 = mod_norm(static_cast<int64_t>(mod_norm(
                                     static_cast<int64_t>(mod_mul(a, rg.sxy2)) -
                                     static_cast<int64_t>(mod_mul(a, rg.sxy1))) -
                                     static_cast<int64_t>(mod_mul(2, rg.sprod2))));
                uint64_t q0 = mod_norm(static_cast<int64_t>(mod_mul(a2,
                                     mod_norm(static_cast<int64_t>(rg.c1) - static_cast<int64_t>(rg.c2)))) +
                                     static_cast<int64_t>(mod_mul(a, rg.sxy2)) -
                                     static_cast<int64_t>(rg.sprod2));

                uint64_t p2 = mod_mul(2, q2);
                uint64_t p1 = mod_mul(2, q1);
                uint64_t p0 = mod_norm(static_cast<int64_t>(mod_mul(2, q0)) +
                                       static_cast<int64_t>(mod_mul(2, a)));

                long long L1 = rg.L1;
                long long L2 = rg.L2;
                uint64_t pow_L1 = pow2_mod(static_cast<unsigned __int128>(L1));
                uint64_t pow_L2p1 = pow2_mod(static_cast<unsigned __int128>(L2) + 1);
                uint64_t s0_r = mod_norm(static_cast<int64_t>(pow_L2p1) - 1);
                uint64_t s0_l = mod_norm(static_cast<int64_t>(pow_L1) - 1);
                uint64_t s0 = mod_norm(static_cast<int64_t>(s0_r) - static_cast<int64_t>(s0_l));

                uint64_t s1_r = sum_lpow_upto(L2, pow_L2p1);
                uint64_t s1_l = sum_lpow_upto(L1 - 1, pow_L1);
                uint64_t s1 = mod_norm(static_cast<int64_t>(s1_r) - static_cast<int64_t>(s1_l));

                uint64_t s2_r = sum_l2pow_upto(L2, pow_L2p1);
                uint64_t s2_l = sum_l2pow_upto(L1 - 1, pow_L1);
                uint64_t s2 = mod_norm(static_cast<int64_t>(s2_r) - static_cast<int64_t>(s2_l));

                uint64_t t1 = sum_range_1(L1, L2);
                uint64_t t2 = sum_range_2(L1, L2);
                uint64_t t3 = sum_range_3(L1, L2);
                uint64_t cnt = static_cast<uint64_t>((L2 - L1 + 1) % MOD);

                uint64_t sum_pow = 0;
                sum_pow = (sum_pow + mod_mul(p2, s2)) % MOD;
                sum_pow = (sum_pow + mod_mul(p1, s1)) % MOD;
                sum_pow = (sum_pow + mod_mul(p0, s0)) % MOD;

                uint64_t sum_poly = 0;
                uint64_t t3_plus_t2 = (t3 + t2) % MOD;
                uint64_t t2_plus_t1 = (t2 + t1) % MOD;
                uint64_t t1_plus_cnt = (t1 + cnt) % MOD;
                sum_poly = (sum_poly + mod_mul(p2, t3_plus_t2)) % MOD;
                sum_poly = (sum_poly + mod_mul(p1, t2_plus_t1)) % MOD;
                sum_poly = (sum_poly + mod_mul(p0, t1_plus_cnt)) % MOD;

                uint64_t term = mod_norm(static_cast<int64_t>(sum_pow) - static_cast<int64_t>(sum_poly));

                local += term;
                if (local >= MOD) local %= MOD;
            }
        }
        partial[idx] = local % MOD;
    };

    if (threads == 1) {
        worker(0);
    } else {
        vector<thread> pool;
        pool.reserve(threads);
        for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();
    }

    uint64_t sum_line = 0;
    for (uint64_t v : partial) {
        sum_line += v;
        if (sum_line >= MOD) sum_line %= MOD;
    }

    uint64_t non_titanic = (1 + points_mod + sum_line) % MOD;
    return mod_norm(static_cast<int64_t>(total_subsets) - static_cast<int64_t>(non_titanic));
}

bool validate(TotientSummatory& ts) {
    struct Test {
        long long n;
        uint64_t expected;
    };
    const Test tests[] = {
        {1, 11},
        {2, 494},
        {4, 33554178},
        {111, 13500401},
        {100000, 63259062},
    };

    for (const auto& test : tests) {
        uint64_t got = compute_T(test.n, ts, 1);
        if (got != test.expected) {
            cerr << "Validation failed for N=" << test.n
                 << ": got " << got << ", expected " << test.expected << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    long long target = DEFAULT_N;
    if (argc > 1) {
        target = max(0LL, atoll(argv[1]));
    }
    unsigned threads = thread::hardware_concurrency();
    if (argc > 2) {
        threads = max(1, atoi(argv[2]));
    }

    long long max_n = max(target, 100000LL);
    TotientSummatory ts(max_n);
    if (!validate(ts)) return 1;

    uint64_t result = compute_T(target, ts, threads);
    cout << result << '\n';
    return 0;
}
