#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
using namespace std;

static constexpr int MOD = 1'000'000'007;

static inline int addmod(int a, int b) { a += b; if (a >= MOD) a -= MOD; return a; }
static inline int submod(int a, int b) { a -= b; if (a < 0) a += MOD; return a; }
static inline int mulmod(long long a, long long b) { return int((a * b) % MOD); }

static int mod_pow(long long a, long long e) {
    long long r = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return (int)r;
}
static int mod_inv(int a) { return mod_pow(a, MOD - 2); }

// inv[i] for 1..n (n < MOD)
static vector<int> compute_inverses(int n) {
    vector<int> inv(n + 1, 0);
    if (n >= 1) inv[1] = 1;
    for (int i = 2; i <= n; i++) {
        inv[i] = int(MOD - (long long)(MOD / i) * inv[MOD % i] % MOD);
    }
    return inv;
}

struct Harmonics {
    vector<int> H;   // H[k] = sum_{i=1..k} 1/i
    vector<int> H2;  // H2[k] = sum_{i=1..k} 1/i^2
};

static Harmonics compute_harmonics(int n, const vector<int>& inv) {
    Harmonics out;
    out.H.assign(n + 1, 0);
    out.H2.assign(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        out.H[i]  = addmod(out.H[i - 1], inv[i]);
        out.H2[i] = addmod(out.H2[i - 1], mulmod(inv[i], inv[i]));
    }
    return out;
}

// Euler phi via linear sieve
static vector<int> compute_phi(int n) {
    vector<int> phi(n + 1, 0);
    vector<int> primes;
    primes.reserve(n / 10);
    vector<bool> isComp(n + 1, false);

    phi[1] = 1;
    for (int i = 2; i <= n; i++) {
        if (!isComp[i]) {
            primes.push_back(i);
            phi[i] = i - 1;
        }
        for (int p : primes) {
            long long v = 1LL * p * i;
            if (v > n) break;
            isComp[(int)v] = true;
            if (i % p == 0) {
                phi[(int)v] = phi[i] * p;
                break;
            } else {
                phi[(int)v] = phi[i] * (p - 1);
            }
        }
    }
    return phi;
}

// S(n) = sum_{a=1..n-1} sum_{b=1..n-a} 1/lcm(a,b)   (mod MOD)
// Uses identity:
//   S(n) = sum_{d=1..n} phi(d)/d^2 * ( H_{floor(n/d)}^2 - H2_{floor(n/d)} )
static int compute_S(int n,
                     const vector<int>& inv,
                     const vector<int>& H,
                     const vector<int>& H2,
                     const vector<int>& phi)
{
    vector<int> pref(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        int inv_i2 = mulmod(inv[i], inv[i]);
        int term = mulmod(phi[i], inv_i2); // phi(i)/i^2
        pref[i] = addmod(pref[i - 1], term);
    }

    auto Tval = [&](int m) -> int {
        long long h = H[m];
        long long t = (h * h) % MOD;
        t = (t - H2[m]) % MOD;
        if (t < 0) t += MOD;
        return (int)t;
    };

    long long S = 0;
    for (int l = 1; l <= n; ) {
        int q = n / l;
        int r = n / q;
        int sumF = submod(pref[r], pref[l - 1]);
        S = (S + 1LL * sumF * Tval(q)) % MOD;
        l = r + 1;
    }
    return (int)S;
}

// Fast solver for given n (requires n < MOD).
static int solve_fast(int n) {
    if (n == 1) return 1;
    if (n == 2) return 5;
    if (n >= MOD) throw runtime_error("n must be < MOD (so denominators are invertible).");

    const auto inv = compute_inverses(n);
    const auto harm = compute_harmonics(n, inv);
    const auto phi = compute_phi(n);
    const int S = compute_S(n, inv, harm.H, harm.H2, phi);

    const int Hn = harm.H[n];
    const int inv2 = inv[2];

    const int inv_n   = inv[n];
    const int inv_nm1 = inv[n - 1];
    const int inv_d   = inv[n - 2]; // n>=3

    // Probabilities (in mod field) for σ = π^i where π uniform in S_n, i uniform in [1..n!]
    const int pF  = mulmod(Hn, inv_n);                                   // P(σ fixes a given point)
    const int pO  = mulmod(submod(1, pF), inv_nm1);                      // P(σ sends x to a specific y!=x)
    const int pSW = mulmod(mulmod(mulmod(harm.H[n / 2], inv2), inv_n), inv_nm1); // P(σ swaps two given points)
    const int pFF = mulmod(mulmod(addmod(submod(n % MOD, Hn), S), inv_n), inv_nm1); // P(σ fixes two given points)

#ifndef NDEBUG
    // Check: for fixed x, P(σ(x)=x) + sum_{y!=x} P(σ(x)=y) == 1
    long long chk = (pF + 1LL * (n - 1) % MOD * pO) % MOD;
    assert(chk == 1);
#endif

    // f(d) = P( σ(a) < σ(a+d) ) is affine in d:
    // f(d) = K0 + K1*d.
    const int K1 = mulmod(addmod(submod(submod(pF, pFF), pO), pSW), inv_d);

    int K0 = pFF;
    K0 = addmod(K0, mulmod(submod(pF, pFF), mulmod((n - 3) % MOD, inv_d)));
    K0 = addmod(K0, mulmod(submod(pO, pSW), mulmod((n - 1) % MOD, inv_d)));
    K0 = addmod(K0, inv2);
    K0 = submod(K0, pF);
    K0 = submod(K0, pO);
    K0 = addmod(K0, mulmod(addmod(pFF, pSW), inv2)); // +(pFF+pSW)/2

    // Factorials
    vector<int> fact(n + 1, 1);
    for (int i = 1; i <= n; i++) fact[i] = mulmod(fact[i - 1], i);

    // Expected Lehmer digit for position j: let m=j-1
    // E[a_{j}] = (n - Hn)/2 + m*((Hn-1)/(n-1) - K0) - K1*m(m+1)/2
    const int term1_const = mulmod(submod(n % MOD, Hn), inv2);          // (n - Hn)/2
    const int A_const     = mulmod(submod(Hn, 1), inv_nm1);             // (Hn - 1)/(n-1)
    const int C1          = submod(A_const, K0);

    // Parallel sum: Sum_{j=1..n} E[a_j] * (n-j)!  where (n-j)! = fact[n-j]
    unsigned T = max(1u, thread::hardware_concurrency());
    if (n < 100'000) T = 1; // avoid overhead for small n

    vector<long long> partial(T, 0);
    vector<thread> pool;
    pool.reserve(T);

    const int block = (n + (int)T - 1) / (int)T;

    for (unsigned t = 0; t < T; t++) {
        int L = (int)t * block + 1;
        int R = min(n, (int)(t + 1) * block);
        if (L > R) continue;

        pool.emplace_back([&, t, L, R]() {
            long long sum = 0;
            for (int j = L; j <= R; j++) {
                long long m = (long long)j - 1;
                long long mm = m % MOD;
                long long tri = mm * ((mm + 1) % MOD) % MOD;
                tri = tri * inv2 % MOD; // m(m+1)/2

                long long Ea = term1_const;
                Ea = (Ea + mm * C1) % MOD;
                Ea = (Ea - 1LL * K1 * tri) % MOD;
                if (Ea < 0) Ea += MOD;

                sum += Ea * fact[n - j] % MOD;
                if (sum > (1LL << 62)) sum %= MOD; // overflow guard
            }
            partial[t] = sum % MOD;
        });
    }

    for (auto &th : pool) th.join();

    long long sum = 0;
    for (auto v : partial) sum = (sum + v) % MOD;

    const int E_rank = addmod(1, (int)sum);
    const long long ans = 1LL * fact[n] * fact[n] % MOD * E_rank % MOD;
    return (int)ans;
}

static int brute_Q_mod(int n) {
    // Only intended for small n (<= 8).
    vector<int> p(n);
    iota(p.begin(), p.end(), 1);
    vector<vector<int>> perms;
    do perms.push_back(p);
    while (next_permutation(p.begin(), p.end()));

    int fact = 1;
    for (int i = 2; i <= n; i++) fact *= i;

    unordered_map<long long, int> rank;
    rank.reserve(perms.size() * 2);

    auto key = [&](const vector<int>& perm) -> long long {
        // n<=8 => digits are single-digit; base-10 concatenation is unique.
        long long k = 0;
        for (int x : perm) k = k * 10 + x;
        return k;
    };
    for (int i = 0; i < (int)perms.size(); i++) rank[key(perms[i])] = i + 1;

    auto compose = [&](const vector<int>& a, const vector<int>& b) -> vector<int> {
        vector<int> c(n);
        for (int i = 0; i < n; i++) c[i] = a[b[i] - 1];
        return c;
    };
    auto lcm_int = [&](int a, int b) -> int { return a / std::gcd(a, b) * b; };
    auto order = [&](const vector<int>& perm) -> int {
        vector<char> vis(n, 0);
        int ord = 1;
        for (int i = 0; i < n; i++) if (!vis[i]) {
            int len = 0, j = i;
            while (!vis[j]) {
                vis[j] = 1;
                j = perm[j] - 1;
                len++;
            }
            ord = lcm_int(ord, len);
        }
        return ord;
    };

    long long total = 0;
    for (const auto &perm : perms) {
        int m = order(perm);
        vector<int> cur = perm;
        long long s = 0;
        for (int k = 1; k <= m; k++) {
            s += rank[key(cur)];
            cur = compose(perm, cur);
        }
        total += 1LL * (fact / m) * s;
    }
    return (int)(total % MOD);
}

static void self_test() {
    // Given checkpoints
    assert(solve_fast(2)  == 5);
    assert(solve_fast(3)  == 88);
    assert(solve_fast(6)  == 133103808);
    assert(solve_fast(10) == 468421536);

    // Brute vs fast for small n
    for (int n = 1; n <= 8; n++) {
        int fast = solve_fast(n);
        int brute = brute_Q_mod(n);
        assert(fast == brute);
    }
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc > 1 && string(argv[1]) == "--test") {
        self_test();
        return 0;
    }

#ifndef NDEBUG
    // Lightweight correctness checkpoints in debug builds
    self_test();
#endif

    int n = 1'000'000;
    if (argc > 1) n = stoi(argv[1]);

    cout << solve_fast(n) << "\n";
    return 0;
}
