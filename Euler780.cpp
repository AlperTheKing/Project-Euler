#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr i64 MOD = 1'000'000'007LL;

static i128 abs_i128(i128 x) {
    return x < 0 ? -x : x;
}

static i128 gcd_i128(i128 a, i128 b) {
    a = abs_i128(a);
    b = abs_i128(b);
    while (b != 0) {
        const i128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static i128 floor_div_i128(const i128 a, const i128 b) {
    i128 q = a / b;
    i128 r = a % b;
    if (r < 0) {
        q -= 1;
    }
    return q;
}

static i128 isqrt_i128(const i128 x) {
    if (x <= 0) {
        return 0;
    }

    long double approx = std::sqrt(static_cast<long double>(x));
    i128 r = static_cast<i128>(approx);

    while ((r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

static i128 floor_sqrt3_mul(const i128 n) {
    return isqrt_i128(3 * n * n);
}

static i128 floor_div_sqrt3(const i128 n, const i128 m) {
    const i128 nn = n * n;
    const i128 den = 3 * m * m;
    i128 t = isqrt_i128(nn / den);
    while (den * (t + 1) * (t + 1) <= nn) {
        ++t;
    }
    while (den * t * t > nn) {
        --t;
    }
    return t;
}

static i128 divisor_summatory(const i64 n) {
    i128 res = 0;
    i64 i = 1;
    while (i <= n) {
        const i64 q = n / i;
        const i64 j = n / q;
        res += static_cast<i128>(q) * static_cast<i128>(j - i + 1);
        i = j + 1;
    }
    return res;
}

struct Alpha {
    i128 a;
    i128 b;
    i128 c;
};

static Alpha normalize_alpha(i128 a, i128 b, i128 c) {
    if (c < 0) {
        a = -a;
        b = -b;
        c = -c;
    }
    const i128 g = gcd_i128(gcd_i128(abs_i128(a), abs_i128(b)), c);
    if (g > 1) {
        a /= g;
        b /= g;
        c /= g;
    }
    return {a, b, c};
}

static i128 floor_qsqrt3(const i128 a, const i128 b, const i128 c) {
    if (b == 0) {
        return floor_div_i128(a, c);
    }

    const i128 fb = floor_sqrt3_mul(b);
    i128 x = floor_div_i128(a + fb, c);
    const i128 bb3 = 3 * b * b;

    while (true) {
        const i128 y = (x + 1) * c - a;
        if (y <= 0 || y * y <= bb3) {
            ++x;
        } else {
            break;
        }
    }

    while (true) {
        const i128 y = x * c - a;
        if (y <= 0 || y * y <= bb3) {
            break;
        }
        --x;
    }

    return x;
}

static i128 alpha_floor(const Alpha& alpha) {
    return floor_qsqrt3(alpha.a, alpha.b, alpha.c);
}

static i128 alpha_mul_floor(const Alpha& alpha, const i128 n) {
    return floor_qsqrt3(alpha.a * n, alpha.b * n, alpha.c);
}

static Alpha alpha_sub_int(const Alpha& alpha, const i128 k) {
    return normalize_alpha(alpha.a - k * alpha.c, alpha.b, alpha.c);
}

static Alpha alpha_div_alpha_minus_one(const Alpha& alpha) {
    const i128 ac = alpha.a - alpha.c;

    i128 A = alpha.a * ac - 3 * alpha.b * alpha.b;
    i128 B = -alpha.b * alpha.c;
    i128 D = ac * ac - 3 * alpha.b * alpha.b;

    if (D < 0) {
        A = -A;
        B = -B;
        D = -D;
    }
    if (B < 0) {
        A = -A;
        B = -B;
    }

    return normalize_alpha(A, B, D);
}

static i128 tri(const i128 n) {
    return n * (n + 1) / 2;
}

static i128 beatty_sum(Alpha alpha, i128 n) {
    i128 res = 0;
    i64 sign = 1;

    while (n > 0) {
        const i128 f = alpha_floor(alpha);
        if (f > 1) {
            res += static_cast<i128>(sign) * (f - 1) * tri(n);
            alpha = alpha_sub_int(alpha, f - 1);
        }

        const i128 m = alpha_mul_floor(alpha, n);
        res += static_cast<i128>(sign) * tri(m);

        n = m - n;
        if (n <= 0) {
            break;
        }

        alpha = alpha_div_alpha_minus_one(alpha);
        sign = -sign;
    }

    return res;
}

static i128 beatty_sqrt3(const i128 c, const i128 n) {
    if (n <= 0) {
        return 0;
    }
    return beatty_sum({0, c, 1}, n);
}

static void linear_sieve(const int n,
                         std::vector<int>& spf,
                         std::vector<int>& mu) {
    spf.assign(n + 1, 0);
    mu.assign(n + 1, 0);

    std::vector<int> primes;
    mu[1] = 1;
    spf[1] = 1;

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
            mu[i] = -1;
        }
        for (const int p : primes) {
            const i64 v = static_cast<i64>(i) * static_cast<i64>(p);
            if (v > n) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
            if (i % p == 0) {
                mu[static_cast<std::size_t>(v)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(v)] = -mu[i];
        }
    }
}

static std::vector<int> distinct_prime_factors(int x,
                                               const std::vector<int>& spf) {
    std::vector<int> ps;
    while (x > 1) {
        const int p = spf[static_cast<std::size_t>(x)];
        ps.push_back(p);
        while (x % p == 0) {
            x /= p;
        }
    }
    return ps;
}

static std::vector<std::vector<std::pair<int, int>>> build_squarefree_divs(
    const int max_n,
    const std::vector<int>& spf
) {
    std::vector<std::vector<std::pair<int, int>>> sf(max_n + 1);
    sf[1] = {{1, 1}};

    for (int n = 2; n <= max_n; ++n) {
        const std::vector<int> ps = distinct_prime_factors(n, spf);
        std::vector<std::pair<int, int>> divs;
        divs.push_back({1, 1});

        for (const int p : ps) {
            const std::size_t cur = divs.size();
            for (std::size_t i = 0; i < cur; ++i) {
                divs.push_back({divs[i].first * p, -divs[i].second});
            }
        }

        sf[n] = std::move(divs);
    }

    return sf;
}

static std::vector<std::vector<int>> build_all_divisors(
    const int max_n,
    const std::vector<int>& spf
) {
    std::vector<std::vector<int>> divs(max_n + 1);
    divs[1] = {1};

    for (int n = 2; n <= max_n; ++n) {
        int x = n;
        const int p = spf[static_cast<std::size_t>(x)];
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }

        const std::vector<int>& base = divs[static_cast<std::size_t>(x)];
        std::vector<int> out;
        out.reserve(base.size() * static_cast<std::size_t>(e + 1));

        int pe = 1;
        for (int i = 0; i <= e; ++i) {
            for (const int d : base) {
                out.push_back(d * pe);
            }
            pe *= p;
        }

        divs[static_cast<std::size_t>(n)] = std::move(out);
    }

    return divs;
}

static unsigned choose_threads(const i64 work_items) {
    if (work_items < 2'000) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 1;
    }

    if (static_cast<i64>(threads) > work_items) {
        threads = static_cast<unsigned>(work_items);
    }
    return std::max(1U, threads);
}

struct StripPartial {
    i128 s1 = 0;
    i128 s2 = 0;
    i128 diag1 = 0;
    i128 diag2 = 0;
};

static std::pair<i128, i128> strip_hyperbola_sum(
    const i64 N,
    const i64 V,
    const i64 L,
    const std::vector<std::vector<std::pair<int, int>>>& sf_divs,
    const std::vector<std::vector<int>>& all_divs
) {
    const unsigned threads = choose_threads(V);
    std::vector<StripPartial> partial(threads);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    std::atomic<i64> next_v{1};

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            StripPartial local;

            while (true) {
                const i64 v = next_v.fetch_add(1, std::memory_order_relaxed);
                if (v > V) {
                    break;
                }

                const i64 Umax = L / v;
                const std::vector<int>& dv = all_divs[static_cast<std::size_t>(v)];

                for (const i64 d : dv) {
                    const i64 m = v / d;
                    const i64 hi = Umax / d;
                    if (hi < m) {
                        continue;
                    }

                    const i64 w = N / (2 * d);
                    const i64 lo_minus = m - 1;

                    i128 cnt = 0;
                    i128 sm = 0;

                    const auto& sfm = sf_divs[static_cast<std::size_t>(m)];
                    for (const auto& [q, muq] : sfm) {
                        cnt += static_cast<i128>(muq) *
                               (static_cast<i128>(hi / q) - static_cast<i128>(lo_minus / q));

                        const i64 hiq = hi / q;
                        const i64 loq = lo_minus / q;
                        if (hiq > 0) {
                            const i128 c = static_cast<i128>(v) * static_cast<i128>(q);
                            sm += static_cast<i128>(muq) *
                                  (beatty_sqrt3(c, hiq) - beatty_sqrt3(c, loq));
                        }
                    }

                    local.s1 += static_cast<i128>(w) * cnt;
                    local.s2 += sm;
                }

                local.diag1 += static_cast<i128>(N / (2 * v));
                local.diag2 += floor_sqrt3_mul(v);
            }

            partial[tid] = local;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    i128 s1 = 0;
    i128 s2 = 0;
    i128 diag1 = 0;
    i128 diag2 = 0;

    for (const StripPartial& p : partial) {
        s1 += p.s1;
        s2 += p.s2;
        diag1 += p.diag1;
        diag2 += p.diag2;
    }

    const i128 S1 = 2 * s1 - diag1;
    const i128 S2 = 2 * s2 - diag2;
    return {S1, S2};
}

static i64 count_mod3_res(const i64 lo, const i64 hi, const i64 r) {
    if (hi < lo) {
        return 0;
    }
    const i64 rem = ((lo % 3) + 3) % 3;
    const i64 delta = (r - rem + 3) % 3;
    const i64 first = lo + delta;
    if (first > hi) {
        return 0;
    }
    return (hi - first) / 3 + 1;
}

static i128 hex_Hsum(
    const i64 X,
    const std::vector<std::vector<std::pair<int, int>>>& sf_divs
) {
    if (X <= 0) {
        return 0;
    }

    std::unordered_map<i64, i128> d_cache;
    d_cache.reserve(1 << 15);

    auto D = [&](const i64 n) -> i128 {
        if (n <= 0) {
            return 0;
        }
        const auto it = d_cache.find(n);
        if (it != d_cache.end()) {
            return it->second;
        }
        const i128 val = divisor_summatory(n);
        d_cache.emplace(n, val);
        return val;
    };

    const i64 V = static_cast<i64>(isqrt_i128(X));
    i128 extra = 0;

    for (i64 v = 1; v <= V; ++v) {
        const i128 disc0 = 4 * static_cast<i128>(X) - 3 * static_cast<i128>(v) * static_cast<i128>(v);
        if (disc0 <= 0) {
            break;
        }

        const i64 Umax = static_cast<i64>((-static_cast<i128>(v) + isqrt_i128(disc0)) / 2);
        i64 u = v + 1;

        const auto& sfv = sf_divs[static_cast<std::size_t>(v)];
        const i64 vmod = v % 3;

        while (u <= Umax) {
            const i128 t = static_cast<i128>(u) * u + static_cast<i128>(u) * v + static_cast<i128>(v) * v;
            const i64 q = static_cast<i64>(X / t);
            if (q == 0) {
                break;
            }

            const i64 T = X / q;
            const i128 disc = 4 * static_cast<i128>(T) - 3 * static_cast<i128>(v) * static_cast<i128>(v);
            i64 uhi = static_cast<i64>((-static_cast<i128>(v) + isqrt_i128(disc)) / 2);
            if (uhi > Umax) {
                uhi = Umax;
            }

            const i64 lo1 = u - 1;

            i64 total = 0;
            for (const auto& [d, mud] : sfv) {
                total += mud * (uhi / d - lo1 / d);
            }

            if (vmod != 0) {
                i64 bad = 0;
                for (const auto& [d, mud] : sfv) {
                    const i64 dm3 = d % 3;
                    const i64 inv = (dm3 == 1 ? 1 : 2);
                    const i64 tlo = (u + d - 1) / d;
                    const i64 thi = uhi / d;
                    const i64 rr = (vmod * inv) % 3;
                    bad += mud * count_mod3_res(tlo, thi, rr);
                }
                total -= bad;
            }

            if (total != 0) {
                extra += static_cast<i128>(total) * D(q);
            }

            u = uhi + 1;
        }
    }

    return D(X) + 2 * extra;
}

static i64 solve_fast(const i64 N) {
    if (N <= 0) {
        return 0;
    }

    const i64 M = N / 2;
    const i64 L = static_cast<i64>(floor_div_sqrt3(M, 1));
    const i64 V_strip = static_cast<i64>(isqrt_i128(L));

    const i64 X = N / 4;
    const i64 V_hex = (X > 0 ? static_cast<i64>(isqrt_i128(X)) : 0);

    const int max_pre = static_cast<int>(std::max<i64>({1, V_strip, V_hex}));

    std::vector<int> spf;
    std::vector<int> mu;
    linear_sieve(max_pre, spf, mu);

    const auto sf_divs = build_squarefree_divs(max_pre, spf);
    const auto all_divs = build_all_divisors(max_pre, spf);

    const i128 base = 2 * divisor_summatory(M);
    const auto [S1, S2] = strip_hyperbola_sum(N, V_strip, L, sf_divs, all_divs);
    const i128 strip_part = base + 4 * (S1 - S2);

    const i128 H = hex_Hsum(X, sf_divs);
    const i128 total = strip_part - 4 * H;

    i128 modded = total % MOD;
    if (modded < 0) {
        modded += MOD;
    }
    return static_cast<i64>(modded);
}

int main() {
    const i64 g5 = solve_fast(5);
    const i64 g6 = solve_fast(6);
    assert(((g6 - g5) % MOD + MOD) % MOD == 8);
    assert(g6 == 14);
    assert(solve_fast(100) == 8'090);
    assert(solve_fast(100'000) == 645'124'048);

    std::cout << solve_fast(1'000'000'000) << '\n';
    return 0;
}
