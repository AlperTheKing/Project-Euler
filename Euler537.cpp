#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using namespace std;

static constexpr uint32_t MOD = 1004535809u; // 2^21 * 479 + 1
static constexpr uint32_t PRIMITIVE_ROOT = 3u;

static uint32_t mod_pow(uint32_t a, uint32_t e) {
    uint64_t r = 1;
    uint64_t x = a;
    while (e) {
        if (e & 1u) r = (r * x) % MOD;
        x = (x * x) % MOD;
        e >>= 1u;
    }
    return (uint32_t)r;
}

static uint32_t mod_inv(uint32_t a) {
    // MOD is prime.
    return mod_pow(a, MOD - 2);
}

static void ntt(vector<uint32_t> &a, bool invert) {
    int n = (int)a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        uint32_t wlen = mod_pow(PRIMITIVE_ROOT, (MOD - 1) / (uint32_t)len);
        if (invert) wlen = mod_inv(wlen);
        for (int i = 0; i < n; i += len) {
            uint64_t w = 1;
            for (int j = 0; j < len / 2; j++) {
                uint32_t u = a[i + j];
                uint32_t v = (uint32_t)((w * a[i + j + len / 2]) % MOD);
                uint32_t x = u + v;
                if (x >= MOD) x -= MOD;
                a[i + j] = x;

                uint32_t y = (u >= v) ? (u - v) : (u + MOD - v);
                a[i + j + len / 2] = y;

                w = (w * wlen) % MOD;
            }
        }
    }

    if (invert) {
        uint32_t inv_n = mod_inv((uint32_t)n);
        for (uint32_t &x : a) x = (uint32_t)((uint64_t)x * inv_n % MOD);
    }
}

static vector<uint32_t> convolution(const vector<uint32_t> &a, const vector<uint32_t> &b, int max_deg) {
    if (a.empty() || b.empty()) return {};
    int need = (int)a.size() + (int)b.size() - 1;
    need = min(need, max_deg + 1);

    int ntt_n = 1;
    while (ntt_n < (int)a.size() + (int)b.size() - 1) ntt_n <<= 1;

    vector<uint32_t> fa(a.begin(), a.end());
    vector<uint32_t> fb(b.begin(), b.end());
    fa.resize(ntt_n);
    fb.resize(ntt_n);

    ntt(fa, false);
    ntt(fb, false);
    for (int i = 0; i < ntt_n; i++) {
        fa[i] = (uint32_t)((uint64_t)fa[i] * fb[i] % MOD);
    }
    ntt(fa, true);

    fa.resize(need);
    return fa;
}

static vector<uint32_t> poly_pow(vector<uint32_t> base, int exp, int max_deg) {
    vector<uint32_t> res(1, 1);
    while (exp > 0) {
        if (exp & 1) res = convolution(res, base, max_deg);
        exp >>= 1;
        if (exp) base = convolution(base, base, max_deg);
    }
    if ((int)res.size() < max_deg + 1) res.resize(max_deg + 1, 0);
    return res;
}

static vector<int> sieve_primes(int n) {
    vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; (int64_t)i * i <= n; i++) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    vector<int> primes;
    for (int i = 2; i <= n; i++) if (is_prime[i]) primes.push_back(i);
    return primes;
}

static vector<uint32_t> build_a(int max_n) {
    // a[0]=1, and for n>=1: a[n]=p_{n+1}-p_n.
    int need_primes = max_n + 1; // p_1..p_{max_n+1}
    int limit = 300000;
    vector<int> primes;
    while (true) {
        primes = sieve_primes(limit);
        if ((int)primes.size() >= need_primes) break;
        limit *= 2;
    }

    vector<uint32_t> a(max_n + 1);
    a[0] = 1;
    for (int n = 1; n <= max_n; n++) {
        a[n] = (uint32_t)(primes[n] - primes[n - 1]);
    }
    return a;
}

static uint32_t solve_T(const vector<uint32_t> &a_full, int n, int k) {
    vector<uint32_t> a(a_full.begin(), a_full.begin() + (n + 1));
    vector<uint32_t> pk = poly_pow(a, k, n);
    return pk[n];
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int N = 20000;
    const int K = 20000;

    vector<uint32_t> a_full = build_a(N);

    // Validation points from the problem statement.
    {
        uint32_t t10 = solve_T(a_full, 10, 10);
        assert(t10 == 869985u);
        uint32_t t1k = solve_T(a_full, 1000, 1000);
        assert(t1k == 578270566u);
    }

    uint32_t ans = solve_T(a_full, N, K);
    cout << ans << "\n";
    return 0;
}
