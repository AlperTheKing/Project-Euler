#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace std;

struct Mod2Pow64 {
    static constexpr int E = 64;
    static constexpr int D = 2 * E - 2;

    array<uint64_t, D + 1> F_small{};

    array<uint64_t, D + 1> dF{};

    array<uint64_t, D + 1> inv_small{};

    static inline uint64_t inv_odd(uint64_t a) {

        if ((a & 1ULL) == 0) {
            throw runtime_error("inv_odd called with even input");
        }
        uint64_t x = 1;
        for (int i = 0; i < 6; ++i) {
            x = x * (uint64_t)(2 - a * x);
        }
        return x;
    }

    Mod2Pow64() {

        F_small[0] = 1;
        for (int u = 1; u <= D; ++u) {
            uint64_t factor = (uint64_t)(2 * (u - 1) + 1);
            F_small[u] = F_small[u - 1] * factor;
        }

        array<uint64_t, D + 1> tmp = F_small;
        for (int k = 0; k <= D; ++k) {
            dF[k] = tmp[0];

            for (int i = 0; i <= D - k - 1; ++i) {
                tmp[i] = tmp[i + 1] - tmp[i];
            }
        }

        inv_small.fill(0);
        for (int x = 1; x <= D; x += 2) {
            inv_small[x] = inv_odd((uint64_t)x);
        }
    }

    inline uint64_t evalF(uint64_t u) const {

        if (u <= (uint64_t)D) return F_small[(size_t)u];

        uint64_t res = dF[0];

        uint64_t binom_odd = 1;
        int exp2 = 0;

        for (int k = 1; k <= D; ++k) {
            uint64_t num = u - (uint64_t)k + 1;
            uint64_t den = (uint64_t)k;

            int tz_num = __builtin_ctzll(num);
            int tz_den = __builtin_ctzll(den);

            exp2 += tz_num - tz_den;

            binom_odd *= (num >> tz_num);
            uint64_t den_odd = (den >> tz_den);
            binom_odd *= inv_small[(size_t)den_odd];

            uint64_t binom = (exp2 >= 64) ? 0ULL : (binom_odd << exp2);
            res += dF[(size_t)k] * binom;
        }
        return res;
    }

    inline uint64_t odd_fact(uint64_t n) const {

        uint64_t res = 1;
        while (n > 0) {
            uint64_t m = (n + 1) >> 1;
            res *= evalF(m);
            n >>= 1;
        }
        return res;
    }

    inline uint64_t odd_central_binom(uint64_t k) const {

        uint64_t a = odd_fact(2 * k);
        uint64_t b = odd_fact(k);
        uint64_t invb = inv_odd(b);
        return a * invb * invb;
    }
};

static inline uint64_t u_value(uint64_t n, const Mod2Pow64& mod) {

    const uint64_t base = n + 1;

    uint64_t O = mod.odd_central_binom(base);
    uint64_t k = base;

    int pc = __builtin_popcountll(k);
    uint64_t B = (pc >= 64) ? 0ULL : (O << pc);

    uint64_t sum = 0;
    for (int j = 0; j < 64; ++j) {
        uint64_t term = B << j;
        if (j & 1) term = (uint64_t)(0) - term;
        sum += term;

        if (j != 63) {

            uint64_t denom = k + 1;
            int tz = __builtin_ctzll(denom);
            uint64_t denom_odd = denom >> tz;
            uint64_t inv = Mod2Pow64::inv_odd(denom_odd);

            O *= (2 * k + 1);
            O *= inv;

            ++k;
            pc = __builtin_popcountll(k);
            B = (pc >= 64) ? 0ULL : (O << pc);
        }
    }

    uint64_t c = (uint64_t)(-3) * sum;

    if (c == 0) {

        throw runtime_error("Need higher 2-adic precision (c==0 mod 2^64)");
    }

    int tz = __builtin_ctzll(c);
    return base + (uint64_t)tz;
}

static inline unsigned long long U_value(int N, const Mod2Pow64& mod) {
    unsigned long long total = 0;
    for (int n = 1; n <= N; ++n) {
        uint64_t x = (uint64_t)n;
        uint64_t cube = x * x * x;
        total += (unsigned long long)u_value(cube, mod);
    }
    return total;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Mod2Pow64 mod;

    auto check = [&](const string& name, unsigned long long got, unsigned long long expected) {
        if (got != expected) {
            cerr << "Validation failed: " << name << " got " << got
                 << " expected " << expected << "\n";
            exit(1);
        }
    };

    check("u(4)",  (unsigned long long)u_value(4, mod), 7ULL);
    check("u(20)", (unsigned long long)u_value(20, mod), 24ULL);
    check("U(5)",  U_value(5, mod), 241ULL);

    const int N = 10000;
    unsigned int T = thread::hardware_concurrency();
    if (T == 0) T = 4;
    if (T > (unsigned int)N) T = N;

    atomic<unsigned long long> global_sum{0};

    int chunk = (N + (int)T - 1) / (int)T;
    vector<thread> threads;
    threads.reserve(T);

    for (unsigned int t = 0; t < T; ++t) {
        int L = (int)t * chunk + 1;
        int R = min(N, (int)(t + 1) * chunk);
        if (L > R) break;

        threads.emplace_back([&, L, R]() {
            unsigned long long local = 0;
            for (int n = L; n <= R; ++n) {
                uint64_t x = (uint64_t)n;
                uint64_t cube = x * x * x;
                local += (unsigned long long)u_value(cube, mod);
            }
            global_sum.fetch_add(local, memory_order_relaxed);
        });
    }

    for (auto& th : threads) th.join();

    unsigned long long ans = global_sum.load(memory_order_relaxed);

    cout << ans << "\n";
    return 0;
}
