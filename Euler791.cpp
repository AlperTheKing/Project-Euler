#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

static constexpr long long MOD = 433494437LL;

static inline long long isqrt_ll(long long x) {
    if (x <= 0) return 0;
    long long r = (long long)std::sqrt((long double)x);
    while ((r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

static inline long long ceil_sqrt_ll(long long x) {
    if (x <= 0) return 0;
    long long r = isqrt_ll(x);
    if (r * r < x) ++r;
    return r;
}

static inline long long sum_interval(long long L, long long H, long long N0) {
    if ((L & 1LL) == 0) ++L;
    if ((H & 1LL) == 0) --H;
    if (L > H) return 0;

    long long cnt = (H - L) / 2 + 1;
    __int128 t = cnt - 1;
    __int128 sum_i = (__int128)cnt * t / 2;
    __int128 sum_i2 = t * (__int128)cnt * (2 * t + 1) / 6;

    __int128 sum_C = (__int128)cnt * L + 2 * sum_i;
    __int128 sum_C2 = (__int128)cnt * L * L + 4 * (__int128)L * sum_i + 4 * sum_i2;

    __int128 sum_num = (__int128)cnt * N0 + sum_C2 + 2 * sum_C;
    __int128 sum_s = sum_num / 2;
    long long mod = (long long)(sum_s % MOD);
    if (mod < 0) mod += MOD;
    return mod;
}

static long long S_value(long long n) {
    if (n <= 0) return 0;
    long long R = 8 * n;

    long long A_max = isqrt_ll(R + 3) - 2;
    if (A_max < 1) return 0;
    if ((A_max & 1LL) == 0) --A_max;

    long long countA = (A_max - 1) / 2 + 1;
    unsigned int T = thread::hardware_concurrency();
    if (T == 0) T = 4;
    if ((long long)T > countA) T = (unsigned int)countA;
    if (T == 0) T = 1;

    atomic<long long> nextA{1};
    vector<long long> partial(T, 0);
    vector<thread> threads;
    threads.reserve(T);

    for (unsigned int t = 0; t < T; ++t) {
        threads.emplace_back([&, t]() {
            long long local = 0;
            for (;;) {
                long long A = nextA.fetch_add(2, memory_order_relaxed);
                if (A > A_max) break;
                long long A2 = A * A;
                long long D = R - A2 - 4 * A - 1;
                if (D < 0) continue;
                long long B_max = isqrt_ll(D) - 2;
                if (B_max < -1) continue;
                if (B_max > A) B_max = A;

                for (long long B = -1; B <= B_max; B += 2) {
                    long long B2 = B * B;
                    long long L = R - A2 - B2 - 4 * A - 4 * B - 5;
                    if (L < 0) continue;
                    long long C_max = isqrt_ll(L);

                    long long C_low = -C_max;
                    long long low2 = -B - 2;
                    if (low2 > C_low) C_low = low2;
                    long long C_high = C_max;
                    if (B < C_high) C_high = B;
                    if (C_low > C_high) continue;

                    long long Lmin = 11 - A2 - B2;
                    long long N0 = A2 + B2 + 2 * A + 2 * B + 3;

                    if (Lmin <= 0) {
                        long long part = sum_interval(C_low, C_high, N0);
                        local += part;
                        if (local >= MOD) local -= MOD;
                    } else {
                        long long cmin = ceil_sqrt_ll(Lmin);
                        if ((cmin & 1LL) == 0) ++cmin;

                        long long neg_high = C_high < -cmin ? C_high : -cmin;
                        if (C_low <= neg_high) {
                            long long part = sum_interval(C_low, neg_high, N0);
                            local += part;
                            if (local >= MOD) local -= MOD;
                        }

                        long long pos_low = C_low > cmin ? C_low : cmin;
                        if (pos_low <= C_high) {
                            long long part = sum_interval(pos_low, C_high, N0);
                            local += part;
                            if (local >= MOD) local -= MOD;
                        }
                    }
                }
            }
            partial[t] = local;
        });
    }

    for (auto& th : threads) th.join();

    long long ans = 0;
    for (long long v : partial) {
        ans += v;
        if (ans >= MOD) ans -= MOD;
    }
    return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    auto check = [&](const string& name, long long got, long long expected) {
        if (got != expected) {
            cerr << "Validation failed: " << name << " got " << got
                 << " expected " << expected << "\n";
            exit(1);
        }
    };

    check("S(5)", S_value(5), 48);
    check("S(10^3)", S_value(1000), 37048340);

    long long ans = S_value(100000000LL);
    cout << ans % MOD << "\n";
    return 0;
}
