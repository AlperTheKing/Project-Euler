#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

using int64 = long long;
using i128 = __int128_t;

int64 pow_int(int64 base, int exp) {
    int64 res = 1;
    for (int i = 0; i < exp; ++i) res *= base;
    return res;
}

int64 mod_pow(int64 base, long long exp, int64 mod) {
    int64 res = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = static_cast<i128>(res) * base % mod;
        base = static_cast<i128>(base) * base % mod;
        exp >>= 1;
    }
    return res;
}

int64 ext_gcd(int64 a, int64 b, int64 &x, int64 &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int64 x1 = 0, y1 = 0;
    int64 g = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return g;
}

int64 mod_inv(int64 a, int64 mod) {
    int64 x = 0, y = 0;
    int64 g = ext_gcd(a, mod, x, y);
    if (g != 1) return 0;
    x %= mod;
    if (x < 0) x += mod;
    return x;
}

int64 compute_mod_prime(long long n, int p) {
    int64 mod_p3 = pow_int(p, 3);
    if (n == 0) return 1 % mod_p3;

    int max_j = static_cast<int>(min<long long>(n, 3LL * p - 1));
    int64 mod_big = pow_int(p, 5); // p^{3+2} is enough since max_j < 3p.

    vector<vector<int64>> binom(max_j + 1, vector<int64>(max_j + 1, 0));
    binom[0][0] = 1;
    for (int j = 1; j <= max_j; ++j) {
        binom[j][0] = 1;
        binom[j][j] = 1;
        for (int i = 1; i < j; ++i) {
            int64 v = binom[j - 1][i - 1] + binom[j - 1][i];
            if (v >= mod_big) v -= mod_big;
            binom[j][i] = v;
        }
    }

    vector<int64> pow_i(max_j + 1, 0);
    for (int i = 0; i <= max_j; ++i) {
        pow_i[i] = mod_pow(i, n, mod_big);
    }

    vector<int> vp_fact(max_j + 1, 0);
    vector<int64> u_fact(max_j + 1, 0);
    u_fact[0] = 1 % mod_p3;
    for (int j = 1; j <= max_j; ++j) {
        int temp = j;
        int cnt = 0;
        while (temp % p == 0) {
            temp /= p;
            ++cnt;
        }
        vp_fact[j] = vp_fact[j - 1] + cnt;
        u_fact[j] = static_cast<i128>(u_fact[j - 1]) * (temp % mod_p3) % mod_p3;
    }

    int64 inv2 = mod_inv(2 % mod_p3, mod_p3);
    int64 pow2 = mod_pow(2, n, mod_p3);
    int64 n_mod = n % mod_p3;
    int64 fall = 1 % mod_p3;
    int64 sum = 0;

    for (int j = 1; j <= max_j; ++j) {
        int64 factor = n_mod - (j - 1);
        factor %= mod_p3;
        if (factor < 0) factor += mod_p3;
        fall = static_cast<i128>(fall) * factor % mod_p3;
        pow2 = static_cast<i128>(pow2) * inv2 % mod_p3; // now 2^{n-j}

        int a = vp_fact[j];
        int64 mod_ext = mod_p3;
        for (int t = 0; t < a; ++t) mod_ext *= p;

        int64 N = 0;
        for (int i = 0; i <= j; ++i) {
            int64 term = static_cast<i128>(binom[j][i]) * pow_i[i] % mod_big;
            if ((j - i) & 1) {
                N -= term;
                if (N < 0) N += mod_big;
            } else {
                N += term;
                if (N >= mod_big) N -= mod_big;
            }
        }
        if (mod_ext != mod_big) N %= mod_ext;
        if (N < 0) N += mod_ext;
        int64 N_div = N;
        for (int t = 0; t < a; ++t) N_div /= p;

        int64 inv_u = mod_inv(u_fact[j], mod_p3);
        int64 stirling = static_cast<i128>(N_div % mod_p3) * inv_u % mod_p3;

        int64 term = static_cast<i128>(stirling) * fall % mod_p3;
        term = static_cast<i128>(term) * pow2 % mod_p3;
        sum += term;
        sum %= mod_p3;
    }

    return sum % mod_p3;
}

unsigned long long compute_naive(unsigned int n) {
    unsigned __int128 sum = 0;
    unsigned __int128 binom = 1;
    for (unsigned int k = 0; k <= n; ++k) {
        if (k > 0) {
            binom = binom * (n - k + 1) / k;
        }
        unsigned __int128 powk = 1;
        for (unsigned int i = 0; i < n; ++i) powk *= k;
        sum += binom * powk;
    }
    return static_cast<unsigned long long>(sum);
}

int64 combine_crt(const vector<int64> &rem, const vector<int64> &mods) {
    int64 M = 1;
    for (int64 m : mods) M *= m;
    int64 result = 0;
    for (size_t i = 0; i < mods.size(); ++i) {
        int64 Mi = M / mods[i];
        int64 inv = mod_inv(Mi % mods[i], mods[i]);
        int64 term = static_cast<i128>(rem[i]) * Mi % M;
        term = static_cast<i128>(term) * inv % M;
        result += term;
        result %= M;
    }
    return result;
}

int64 compute_mod_all(long long n) {
    vector<int> primes = {83, 89, 97};
    vector<int64> mods(3, 0);
    vector<int64> rems(3, 0);

    vector<thread> threads;
    threads.reserve(primes.size());
    for (size_t i = 0; i < primes.size(); ++i) {
        threads.emplace_back([&, i]() {
            int p = primes[i];
            mods[i] = pow_int(p, 3);
            rems[i] = compute_mod_prime(n, p);
        });
    }
    for (auto &t : threads) t.join();

    return combine_crt(rems, mods);
}

} // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const unsigned long long s10 = compute_naive(10);
    if (s10 != 142469423360ULL) {
        cerr << "Validation failed: S(10) != 142469423360\n";
        return 1;
    }

    int64 mod_all = pow_int(83, 3) * pow_int(89, 3) * pow_int(97, 3);
    int64 s10_mod = static_cast<int64>(s10 % mod_all);
    int64 s10_fast = compute_mod_all(10);
    if (s10_fast != s10_mod) {
        cerr << "Validation failed: S(10) mod M mismatch\n";
        return 1;
    }

    long long n = 1'000'000'000'000'000'000LL;
    int64 ans = compute_mod_all(n);
    cout << ans << "\n";
    return 0;
}
