#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr int MOD = 1000000007;
constexpr uint64_t MAIN_N = 1000000000000ULL;
constexpr int K_MAX = 50;

inline int addmod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

inline int submod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

inline int mulmod(long long a, long long b) {
    return static_cast<int>((static_cast<__int128>(a) * b) % MOD);
}

struct PowerSum {
    int k = 0;
    int m = 0;
    vector<int> y;
    vector<int> invfact;

    explicit PowerSum(int kk) : k(kk) {
        m = k + 2;
        y.assign(m + 1, 0);
        for (int i = 1; i <= m; ++i) {
            long long p = 1;
            long long b = i;
            int e = k;
            while (e > 0) {
                if (e & 1) p = (p * b) % MOD;
                b = (b * b) % MOD;
                e >>= 1;
            }
            y[i] = addmod(y[i - 1], static_cast<int>(p));
        }

        vector<int> fact(m + 1, 1);
        invfact.assign(m + 1, 1);
        for (int i = 1; i <= m; ++i) fact[i] = mulmod(fact[i - 1], i);

        auto modpow = [&](long long a, long long e) {
            long long r = 1 % MOD;
            a %= MOD;
            while (e > 0) {
                if (e & 1) r = (r * a) % MOD;
                a = (a * a) % MOD;
                e >>= 1;
            }
            return static_cast<int>(r);
        };

        invfact[m] = modpow(fact[m], MOD - 2);
        for (int i = m; i >= 1; --i) invfact[i - 1] = mulmod(invfact[i], i);
    }

    int eval(uint64_t n) const {
        if (n <= static_cast<uint64_t>(m)) return y[static_cast<size_t>(n)];

        const int nm = static_cast<int>(n % MOD);
        vector<int> pre(m + 2, 0), suf(m + 2, 0);

        pre[0] = 1;
        for (int i = 0; i <= m; ++i) {
            pre[i + 1] = mulmod(pre[i], (nm - i + MOD) % MOD);
        }

        suf[m + 1] = 1;
        for (int i = m; i >= 0; --i) {
            suf[i] = mulmod(suf[i + 1], (nm - i + MOD) % MOD);
        }

        int ans = 0;
        for (int i = 0; i <= m; ++i) {
            int num = mulmod(pre[i], suf[i + 1]);
            int den = mulmod(invfact[i], invfact[m - i]);
            int term = mulmod(y[i], mulmod(num, den));
            if ((m - i) & 1) ans = submod(ans, term);
            else ans = addmod(ans, term);
        }
        return ans;
    }
};

struct SolverK {
    uint64_t N;
    int k;
    int sq;
    uint32_t next_prime_after_sq;

    vector<int> primes;       // primes <= sqrt(N)
    vector<uint64_t> vals;    // distinct values floor(N / i), descending
    vector<int> id_small;
    vector<int> id_large;
    vector<int> g;            // sum_{p<=vals[i]} p^k
    vector<int> prime_prefix; // prefix sums of p^k over primes <= sqrt(N)
    vector<int> ppk;

    explicit SolverK(uint64_t n, int kk)
        : N(n), k(kk), sq(0), next_prime_after_sq(0) {
        sq = static_cast<int>(sqrt(static_cast<long double>(N)));
        while (static_cast<uint64_t>(sq + 1) * static_cast<uint64_t>(sq + 1) <= N) ++sq;
        while (static_cast<uint64_t>(sq) * static_cast<uint64_t>(sq) > N) --sq;

        sieve_primes();
        build_vals();
        build_prime_sums();
    }

    void sieve_primes() {
        vector<bool> comp(sq + 1, false);
        for (int i = 2; i <= sq; ++i) {
            if (comp[i]) continue;
            primes.push_back(i);
            if (1LL * i * i <= sq) {
                for (long long j = 1LL * i * i; j <= sq; j += i) {
                    comp[static_cast<size_t>(j)] = true;
                }
            }
        }

        auto is_prime_u32 = [&](uint32_t x) {
            if (x < 2) return false;
            if ((x & 1U) == 0U) return x == 2U;
            uint32_t r = static_cast<uint32_t>(sqrt(static_cast<long double>(x)));
            for (uint32_t d = 3; d <= r; d += 2) {
                if (x % d == 0U) return false;
            }
            return true;
        };

        uint32_t cand = static_cast<uint32_t>(sq + 1);
        while (!is_prime_u32(cand)) ++cand;
        next_prime_after_sq = cand;
    }

    void build_vals() {
        vals.reserve(static_cast<size_t>(2.0L * sqrt(static_cast<long double>(N)) + 8));
        id_small.assign(sq + 1, -1);
        id_large.assign(sq + 1, -1);

        for (uint64_t l = 1, r; l <= N; l = r + 1) {
            uint64_t v = N / l;
            r = N / v;
            int idx = static_cast<int>(vals.size());
            vals.push_back(v);
            if (v <= static_cast<uint64_t>(sq)) id_small[static_cast<size_t>(v)] = idx;
            else id_large[static_cast<size_t>(N / v)] = idx;
        }
    }

    inline int id(uint64_t x) const {
        if (x <= static_cast<uint64_t>(sq)) return id_small[static_cast<size_t>(x)];
        return id_large[static_cast<size_t>(N / x)];
    }

    static int modpow_int(long long a, int e) {
        long long r = 1 % MOD;
        a %= MOD;
        while (e > 0) {
            if (e & 1) r = (r * a) % MOD;
            a = (a * a) % MOD;
            e >>= 1;
        }
        return static_cast<int>(r);
    }

    void build_prime_sums() {
        const int m = static_cast<int>(vals.size());
        g.assign(m, 0);

        PowerSum ps(k);
        for (int i = 0; i < m; ++i) {
            g[i] = submod(ps.eval(vals[i]), 1);
        }

        const int np = static_cast<int>(primes.size());
        ppk.assign(np, 0);
        prime_prefix.assign(np + 1, 0);
        for (int i = 0; i < np; ++i) {
            ppk[i] = modpow_int(primes[i], k);
            prime_prefix[i + 1] = addmod(prime_prefix[i], ppk[i]);
        }

        for (int i = 0; i < np; ++i) {
            uint64_t p = static_cast<uint64_t>(primes[i]);
            uint64_t p2 = p * p;
            if (p2 > N) break;

            int pk = ppk[i];
            int pref = prime_prefix[i];
            for (int j = 0; j < m && vals[j] >= p2; ++j) {
                int idx = id(vals[j] / p);
                int delta = submod(g[idx], pref);
                g[j] = submod(g[j], mulmod(pk, delta));
            }
        }
    }

    inline int prime_sum(uint64_t x) const {
        if (x < 2) return 0;
        return g[id(x)];
    }

    int dfs(uint64_t n, int start) {
        if (n < 2) return 0;

        const int np = static_cast<int>(primes.size());
        if (start < np) {
            if (static_cast<uint64_t>(primes[start]) > n) return 0;
        } else {
            if (static_cast<uint64_t>(next_prime_after_sq) > n) return 0;
        }

        int ans = submod(prime_sum(n), prime_prefix[start]);

        for (int i = start; i < np; ++i) {
            uint64_t p = static_cast<uint64_t>(primes[i]);
            if (p * p > n) break;

            int pk = ppk[i];
            uint64_t pe = p;
            int e = 1;
            while (pe <= n) {
                ans = addmod(ans, mulmod(pk, dfs(n / pe, i + 1)));
                if (e >= 2) ans = addmod(ans, pk);
                if (pe > n / p) break;
                pe *= p;
                ++e;
            }
        }
        return ans;
    }

    int solve() {
        return addmod(1, dfs(N, 0));
    }
};

bool run_validations() {
    struct CheckK {
        uint64_t n;
        int k;
        int expected;
    };

    const vector<CheckK> checks_k = {
        {10ULL, 1, 41},
        {100ULL, 1, 3512},
        {100ULL, 2, 208090},
        {10000ULL, 1, 35252550},
    };

    for (const auto& c : checks_k) {
        SolverK s(c.n, c.k);
        int got = s.solve();
        if (got != c.expected) {
            cerr << "Validation failed: S_" << c.k << "(" << c.n << ")="
                 << got << ", expected " << c.expected << '\n';
            return false;
        }
    }

    int sum = 0;
    for (int k = 1; k <= 3; ++k) {
        SolverK s(100000000ULL, k);
        sum = addmod(sum, s.solve());
    }
    if (sum != 338787512) {
        cerr << "Validation failed: sum_{k=1..3} S_k(1e8)="
             << sum << ", expected 338787512\n";
        return false;
    }

    cerr << "All validations passed.\n";
    return true;
}

int solve_main(int threads) {
    if (threads <= 0) {
        unsigned hc = thread::hardware_concurrency();
        threads = (hc == 0U) ? 4 : static_cast<int>(hc);
    }
    threads = max(1, min(threads, K_MAX));

    vector<int> sk(K_MAX + 1, 0);
    atomic<int> next_k(1);

    auto worker = [&]() {
        while (true) {
            int k = next_k.fetch_add(1, memory_order_relaxed);
            if (k > K_MAX) break;
            SolverK solver(MAIN_N, k);
            sk[k] = solver.solve();
        }
    };

    vector<thread> pool;
    pool.reserve(static_cast<size_t>(threads));
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    int ans = 0;
    for (int k = 1; k <= K_MAX; ++k) ans = addmod(ans, sk[k]);
    return ans;
}

}  // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int threads = 0;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg.rfind("--threads=", 0) == 0) {
            threads = stoi(arg.substr(10));
        } else if (arg == "--no-validate") {
            validate = false;
        }
    }

    if (validate && !run_validations()) return 1;

    auto start = chrono::high_resolution_clock::now();
    int ans = solve_main(threads);
    auto end = chrono::high_resolution_clock::now();

    cerr << "Elapsed: " << chrono::duration<double>(end - start).count() << "s\n";
    cout << ans << '\n';
    return 0;
}
