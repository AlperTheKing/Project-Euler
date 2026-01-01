#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr long long kMod = 1'000'000'007LL;

long long mod_add(long long a, long long b) {
    a += b;
    if (a >= kMod) a -= kMod;
    return a;
}

long long mod_sub(long long a, long long b) {
    a -= b;
    if (a < 0) a += kMod;
    return a;
}

long long mod_mul(long long a, long long b) {
    return static_cast<long long>((__int128)a * b % kMod);
}

long long mod_pow(long long a, long long e) {
    long long r = 1 % kMod;
    a %= kMod;
    if (a < 0) a += kMod;
    while (e > 0) {
        if (e & 1) r = mod_mul(r, a);
        a = mod_mul(a, a);
        e >>= 1;
    }
    return r;
}

long long mod_inv(long long a) {
    return mod_pow(a, kMod - 2);
}

vector<int> sieve_primes(int n) {
    vector<bool> is_prime(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int i = 2; i * 1LL * i <= n; ++i) {
        if (is_prime[i]) {
            for (long long j = 1LL * i * i; j <= n; j += i) {
                is_prime[static_cast<size_t>(j)] = false;
            }
        }
    }
    vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

int vp_factorial(int n, int p) {
    int e = 0;
    while (n) {
        n /= p;
        e += n;
    }
    return e;
}

long long factorial_mod(int n) {
    long long r = 1;
    for (int i = 2; i <= n; ++i) {
        r = mod_mul(r, i);
    }
    return r;
}

long long geom_sum(long long base, long long e) {
    base %= kMod;
    if (base < 0) base += kMod;
    if (e < 0) return 0;
    if (base == 1) return (e + 1) % kMod;
    long long num = mod_sub(mod_pow(base, e + 1), 1);
    long long den = mod_sub(base, 1);
    return mod_mul(num, mod_inv(den));
}

long long sigma_k_factorial(const vector<pair<int, int>>& pe, int k) {
    long long res = 1;
    for (auto [p, e] : pe) {
        long long base = mod_pow(p, k);
        long long s = geom_sum(base, e);
        res = mod_mul(res, s);
    }
    return res;
}

struct Mat2 {
    long long a00, a01, a10, a11;
};

Mat2 mat_mul(const Mat2& A, const Mat2& B) {
    Mat2 C;
    C.a00 = (mod_mul(A.a00, B.a00) + mod_mul(A.a01, B.a10)) % kMod;
    C.a01 = (mod_mul(A.a00, B.a01) + mod_mul(A.a01, B.a11)) % kMod;
    C.a10 = (mod_mul(A.a10, B.a00) + mod_mul(A.a11, B.a10)) % kMod;
    C.a11 = (mod_mul(A.a10, B.a01) + mod_mul(A.a11, B.a11)) % kMod;
    return C;
}

Mat2 mat_pow(Mat2 base, long long e) {
    Mat2 r{1, 0, 0, 1};
    while (e > 0) {
        if (e & 1) r = mat_mul(r, base);
        base = mat_mul(base, base);
        e >>= 1;
    }
    return r;
}

long long tau_prime_power(long long tau_p, int p, int e) {
    if (e == 0) return 1;
    if (e == 1) return tau_p;
    long long p11 = mod_pow(p, 11);
    Mat2 M{tau_p, mod_sub(0, p11), 1, 0};
    Mat2 P = mat_pow(M, e - 1);
    return (mod_mul(P.a00, tau_p) + mod_mul(P.a01, 1)) % kMod;
}

vector<long long> compute_sigma1_upto(int N) {
    vector<long long> sig(N + 1, 0);
    for (int d = 1; d <= N; ++d) {
        for (int m = d; m <= N; m += d) {
            sig[m] += d;
        }
    }
    for (int i = 0; i <= N; ++i) sig[i] %= kMod;
    return sig;
}

vector<long long> compute_tau_upto(int N, const vector<long long>& sigma1) {
    vector<long long> tau(N + 1, 0);
    tau[0] = 0;
    tau[1] = 1;
    const long long neg24 = mod_sub(0, 24);
    for (int n = 2; n <= N; ++n) {
        long long s = 0;
        for (int k = 1; k < n; ++k) {
            s = (s + sigma1[k] * tau[n - k]) % kMod;
        }
        long long inv = mod_inv(n - 1);
        tau[n] = mod_mul(mod_mul(neg24, s), inv);
    }
    return tau;
}

long long compute_tau_factorial_parallel(const vector<pair<int, int>>& pe,
                                         const vector<long long>& tau_small) {
    unsigned hw = thread::hardware_concurrency();
    unsigned T = hw ? hw : 4;
    if (T > 8) T = 8;

    vector<long long> partial(T, 1);
    vector<thread> threads;
    threads.reserve(T);
    int m = static_cast<int>(pe.size());

    auto worker = [&](unsigned tid) {
        int L = static_cast<int>((long long)m * tid / T);
        int R = static_cast<int>((long long)m * (tid + 1) / T);
        long long prod = 1;
        for (int idx = L; idx < R; ++idx) {
            int p = pe[idx].first;
            int e = pe[idx].second;
            long long tau_p = tau_small[p];
            long long t = tau_prime_power(tau_p, p, e);
            prod = mod_mul(prod, t);
        }
        partial[tid] = prod;
    };

    for (unsigned t = 0; t < T; ++t) threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    long long total = 1;
    for (unsigned t = 0; t < T; ++t) total = mod_mul(total, partial[t]);
    return total;
}

long long conv_Rn_small(int n, int M) {
    vector<long long> sig1 = compute_sigma1_upto(M);
    vector<long long> f(M + 1, 0);
    for (int m = 1; m <= M; ++m) {
        f[m] = (2 * sig1[m]) % kMod;
    }
    vector<long long> dp(M + 1, 0), ndp(M + 1, 0);
    dp[0] = 1;
    for (int iter = 0; iter < n; ++iter) {
        fill(ndp.begin(), ndp.end(), 0);
        for (int s = 0; s <= M; ++s) {
            if (!dp[s]) continue;
            for (int m = 1; s + m <= M; ++m) {
                if (!f[m]) continue;
                ndp[s + m] = (ndp[s + m] + dp[s] * f[m]) % kMod;
            }
        }
        dp.swap(ndp);
    }
    return dp[M];
}

}  // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int N = 10000;

    vector<int> primes = sieve_primes(N);
    vector<pair<int, int>> pe;
    pe.reserve(primes.size());
    for (int p : primes) {
        pe.push_back({p, vp_factorial(N, p)});
    }

    long long n_mod = factorial_mod(N);
    long long n_pow[6];
    n_pow[0] = 1;
    for (int i = 1; i <= 5; ++i) n_pow[i] = mod_mul(n_pow[i - 1], n_mod);

    map<int, long long> sig;
    vector<int> ks = {1, 3, 5, 7, 9, 11};
    mutex sig_mtx;
    vector<thread> sig_threads;
    for (int k : ks) {
        sig_threads.emplace_back([&, k]() {
            long long v = sigma_k_factorial(pe, k);
            lock_guard<mutex> lk(sig_mtx);
            sig[k] = v;
        });
    }
    for (auto& th : sig_threads) th.join();

    vector<long long> sigma1_small = compute_sigma1_upto(N);
    vector<long long> tau_small = compute_tau_upto(N, sigma1_small);

    auto norm = [](long long x) {
        x %= kMod;
        if (x < 0) x += kMod;
        return x;
    };

    if (norm(tau_small[1]) != 1 || norm(tau_small[2]) != norm(-24) ||
        norm(tau_small[3]) != 252 || norm(tau_small[5]) != 4830) {
        cerr << "[Validation failed] tau initial values mismatch\n";
        return 1;
    }

    long long tau_fact = compute_tau_factorial_parallel(pe, tau_small);

    long long c2 = mod_mul(norm(-24), sig[1]);
    long long c4 = mod_mul(240, sig[3]);
    long long c6 = mod_mul(norm(-504), sig[5]);
    long long c8 = mod_mul(480, sig[7]);
    long long c10 = mod_mul(norm(-264), sig[9]);
    long long factor12 = mod_mul(65520, mod_inv(691));
    long long c12 = mod_mul(factor12, sig[11]);

    long long inv2 = mod_inv(2);
    long long inv5 = mod_inv(5);
    long long inv7 = mod_inv(7);
    long long inv24185 = mod_inv(24185);

    long long a1 = c2;

    long long a2 = mod_add(c4, mod_mul(12, mod_mul(n_pow[1], c2)));

    long long a3 = c6;
    a3 = mod_add(a3, mod_mul(9, mod_mul(n_pow[1], c4)));
    a3 = mod_add(a3, mod_mul(72, mod_mul(n_pow[2], c2)));

    long long a4 = c8;
    a4 = mod_add(a4, mod_mul(8, mod_mul(n_pow[1], c6)));
    a4 = mod_add(a4, mod_mul(mod_mul(216, inv5), mod_mul(n_pow[2], c4)));
    a4 = mod_add(a4, mod_mul(288, mod_mul(n_pow[3], c2)));

    long long a5 = c10;
    a5 = mod_add(a5, mod_mul(mod_mul(15, inv2), mod_mul(n_pow[1], c8)));
    a5 = mod_add(a5, mod_mul(mod_mul(240, inv7), mod_mul(n_pow[2], c6)));
    a5 = mod_add(a5, mod_mul(144, mod_mul(n_pow[3], c4)));
    a5 = mod_add(a5, mod_mul(864, mod_mul(n_pow[4], c2)));

    long long a6 = c12;
    a6 = mod_add(a6, mod_mul(mod_mul(norm(-4608), inv24185), tau_fact));
    a6 = mod_add(a6, mod_mul(mod_mul(36, inv5), mod_mul(n_pow[1], c10)));
    a6 = mod_add(a6, mod_mul(30, mod_mul(n_pow[2], c8)));
    a6 = mod_add(a6, mod_mul(mod_mul(720, inv7), mod_mul(n_pow[3], c6)));
    a6 = mod_add(a6, mod_mul(mod_mul(2592, inv7), mod_mul(n_pow[4], c4)));
    a6 = mod_add(a6, mod_mul(mod_mul(10368, inv5), mod_mul(n_pow[5], c2)));

    long long S = 0;
    S = mod_add(S, mod_mul(norm(-6), a1));
    S = mod_add(S, mod_mul(15, a2));
    S = mod_add(S, mod_mul(norm(-20), a3));
    S = mod_add(S, mod_mul(15, a4));
    S = mod_add(S, mod_mul(norm(-6), a5));
    S = mod_add(S, a6);

    long long inv2985984 = mod_inv(2985984);
    long long ans = mod_mul(S, inv2985984);

    {
        long long r1_10 = (2 * (1 + 2 + 5 + 10)) % kMod;
        if (r1_10 != 36) {
            cerr << "[Validation failed] R1(10) expected 36\n";
            return 1;
        }
        long long r2_100 = conv_Rn_small(2, 100);
        if (r2_100 != 1873044) {
            cerr << "[Validation failed] R2(100) expected 1873044, got " << r2_100
                 << "\n";
            return 1;
        }
        long long r6_20_brut = conv_Rn_small(6, 20);

        auto sigma_small = [&](int n, int k) -> long long {
            long long res = 0;
            for (int d = 1; d * 1LL * d <= n; ++d) {
                if (n % d != 0) continue;
                res = (res + mod_pow(d, k)) % kMod;
                if (d * 1LL * d != n) {
                    res = (res + mod_pow(n / d, k)) % kMod;
                }
            }
            return res;
        };
        auto coeffE = [&](int n, int w) -> long long {
            if (w == 2) return mod_mul(norm(-24), sigma_small(n, 1));
            if (w == 4) return mod_mul(240, sigma_small(n, 3));
            if (w == 6) return mod_mul(norm(-504), sigma_small(n, 5));
            if (w == 8) return mod_mul(480, sigma_small(n, 7));
            if (w == 10) return mod_mul(norm(-264), sigma_small(n, 9));
            if (w == 12) return mod_mul(factor12, sigma_small(n, 11));
            return 0LL;
        };

        long long tau20 = tau_small[20];
        long long n20 = 20;
        long long n20p[6];
        n20p[0] = 1;
        for (int i = 1; i <= 5; ++i) n20p[i] = mod_mul(n20p[i - 1], n20);

        long long c2_20 = coeffE(20, 2);
        long long c4_20 = coeffE(20, 4);
        long long c6_20 = coeffE(20, 6);
        long long c8_20 = coeffE(20, 8);
        long long c10_20 = coeffE(20, 10);
        long long c12_20 = coeffE(20, 12);

        long long a1_20 = c2_20;
        long long a2_20 = mod_add(c4_20, mod_mul(12, mod_mul(n20p[1], c2_20)));

        long long a3_20 = c6_20;
        a3_20 = mod_add(a3_20, mod_mul(9, mod_mul(n20p[1], c4_20)));
        a3_20 = mod_add(a3_20, mod_mul(72, mod_mul(n20p[2], c2_20)));

        long long a4_20 = c8_20;
        a4_20 = mod_add(a4_20, mod_mul(8, mod_mul(n20p[1], c6_20)));
        a4_20 = mod_add(a4_20, mod_mul(mod_mul(216, inv5), mod_mul(n20p[2], c4_20)));
        a4_20 = mod_add(a4_20, mod_mul(288, mod_mul(n20p[3], c2_20)));

        long long a5_20 = c10_20;
        a5_20 =
            mod_add(a5_20, mod_mul(mod_mul(15, inv2), mod_mul(n20p[1], c8_20)));
        a5_20 =
            mod_add(a5_20, mod_mul(mod_mul(240, inv7), mod_mul(n20p[2], c6_20)));
        a5_20 = mod_add(a5_20, mod_mul(144, mod_mul(n20p[3], c4_20)));
        a5_20 = mod_add(a5_20, mod_mul(864, mod_mul(n20p[4], c2_20)));

        long long a6_20 = c12_20;
        a6_20 = mod_add(a6_20, mod_mul(mod_mul(norm(-4608), inv24185), tau20));
        a6_20 =
            mod_add(a6_20, mod_mul(mod_mul(36, inv5), mod_mul(n20p[1], c10_20)));
        a6_20 = mod_add(a6_20, mod_mul(30, mod_mul(n20p[2], c8_20)));
        a6_20 =
            mod_add(a6_20, mod_mul(mod_mul(720, inv7), mod_mul(n20p[3], c6_20)));
        a6_20 =
            mod_add(a6_20, mod_mul(mod_mul(2592, inv7), mod_mul(n20p[4], c4_20)));
        a6_20 =
            mod_add(a6_20, mod_mul(mod_mul(10368, inv5), mod_mul(n20p[5], c2_20)));

        long long S20 = 0;
        S20 = mod_add(S20, mod_mul(norm(-6), a1_20));
        S20 = mod_add(S20, mod_mul(15, a2_20));
        S20 = mod_add(S20, mod_mul(norm(-20), a3_20));
        S20 = mod_add(S20, mod_mul(15, a4_20));
        S20 = mod_add(S20, mod_mul(norm(-6), a5_20));
        S20 = mod_add(S20, a6_20);

        long long r6_20_form = mod_mul(S20, inv2985984);
        if (r6_20_form != r6_20_brut) {
            cerr << "[Validation failed] R6(20) mismatch: brut=" << r6_20_brut
                 << " form=" << r6_20_form << "\n";
            return 1;
        }
    }

    cout << ans << "\n";
    return 0;
}
