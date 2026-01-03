#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;

static constexpr int MOD = 1'000'000'007;
static constexpr int K = 30;

// Count labeled tuples using inclusion-exclusion on equal bases; each partition
// yields block sums S and a generating function 1/∏(1 - x^S). Then divide by
// 2!3!4! to forget permutations of identical exponents.

static inline int addmod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static inline int submod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

static int mod_pow(int base, long long exp) {
    long long res = 1;
    long long b = base % MOD;
    while (exp > 0) {
        if (exp & 1LL) res = (res * b) % MOD;
        b = (b * b) % MOD;
        exp >>= 1LL;
    }
    return static_cast<int>(res);
}

struct Entry {
    vector<int> sums;
    int weight_mod;
    int gcd;
};

static vector<Entry> build_entries() {
    const vector<int> exps = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};
    const int n = static_cast<int>(exps.size());

    long long fact[11];
    fact[0] = 1;
    for (int i = 1; i <= 10; ++i) fact[i] = fact[i - 1] * i;

    map<vector<int>, long long> weight_map;
    vector<int> sums;
    vector<int> sizes;

    auto dfs = [&](auto&& self, int idx) -> void {
        if (idx == n) {
            vector<int> key = sums;
            sort(key.begin(), key.end());
            int blocks = static_cast<int>(sizes.size());
            long long mu = ((n - blocks) % 2 == 0) ? 1 : -1;
            for (int sz : sizes) mu *= fact[sz - 1];
            weight_map[key] += mu;
            return;
        }

        for (size_t i = 0; i < sums.size(); ++i) {
            sums[i] += exps[idx];
            sizes[i] += 1;
            self(self, idx + 1);
            sizes[i] -= 1;
            sums[i] -= exps[idx];
        }

        sums.push_back(exps[idx]);
        sizes.push_back(1);
        self(self, idx + 1);
        sums.pop_back();
        sizes.pop_back();
    };

    dfs(dfs, 0);

    vector<Entry> entries;
    entries.reserve(weight_map.size());
    for (const auto& kv : weight_map) {
        long long w = kv.second;
        if (w == 0) continue;
        int g = 0;
        for (int s : kv.first) g = std::gcd(g, s);
        int wmod = static_cast<int>(w % MOD);
        if (wmod < 0) wmod += MOD;
        entries.push_back({kv.first, wmod, g});
    }
    return entries;
}

static const vector<Entry>& get_entries() {
    static vector<Entry> entries = build_entries();
    return entries;
}

static vector<int> sieve_primes(int n) {
    vector<int> primes;
    if (n < 2) return primes;
    vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i <= n; ++i) {
        if (!is_prime[i]) continue;
        primes.push_back(i);
        if (static_cast<long long>(i) * i <= n) {
            for (int j = i * i; j <= n; j += i) is_prime[j] = false;
        }
    }
    return primes;
}

struct ExpCounts {
    vector<pair<int, int>> small;
    vector<pair<int, int>> big;
    int maxSmall;
    bool hasExp1;
};

static ExpCounts build_exp_counts(int n) {
    vector<int> primes = sieve_primes(n);
    vector<int> cnt(n + 1, 0);
    for (int p : primes) {
        int e = 0;
        int m = n;
        while (m > 0) {
            m /= p;
            e += m;
        }
        cnt[e] += 1;
    }

    int maxSmall = max(30, static_cast<int>(sqrt(static_cast<double>(n))));

    vector<pair<int, int>> small;
    vector<pair<int, int>> big;
    for (int e = 1; e <= n; ++e) {
        if (cnt[e] == 0) continue;
        if (e <= maxSmall) {
            small.emplace_back(e, cnt[e]);
        } else {
            big.emplace_back(e, cnt[e]);
        }
    }

    return {small, big, maxSmall, (n >= 2 ? cnt[1] > 0 : false)};
}

using Vec = array<int, K>;

static Vec combine_vecs(const Vec& a, const Vec& b, const Vec& rec) {
    long long tmp[2 * K] = {};
    for (int i = 0; i < K; ++i) {
        if (a[i] == 0) continue;
        for (int j = 0; j < K; ++j) {
            if (b[j] == 0) continue;
            tmp[i + j] = (tmp[i + j] + 1LL * a[i] * b[j]) % MOD;
        }
    }

    for (int i = 2 * K - 2; i >= K; --i) {
        long long val = tmp[i];
        if (val == 0) continue;
        for (int j = 1; j <= K; ++j) {
            tmp[i - j] = (tmp[i - j] + val * rec[j - 1]) % MOD;
        }
    }

    Vec res{};
    for (int i = 0; i < K; ++i) res[i] = static_cast<int>(tmp[i]);
    return res;
}

static int linear_rec(const Vec& rec, const int* init, long long n) {
    if (n < K) return init[n];
    Vec pol{};
    Vec base{};
    pol[0] = 1;
    base[1] = 1; // x

    long long m = n;
    while (m > 0) {
        if (m & 1LL) pol = combine_vecs(pol, base, rec);
        base = combine_vecs(base, base, rec);
        m >>= 1LL;
    }

    long long ans = 0;
    for (int i = 0; i < K; ++i) {
        ans = (ans + 1LL * pol[i] * init[i]) % MOD;
    }
    return static_cast<int>(ans);
}

static int compute_contribution(const Entry& entry,
                                const vector<pair<int, int>>& small,
                                const vector<pair<int, int>>& big,
                                int maxSmall,
                                vector<int>& dp) {
    // Unbounded knapsack for g(e) up to sqrt(n); bigger exponents use the recurrence.
    fill(dp.begin(), dp.end(), 0);
    dp[0] = 1;
    for (int s : entry.sums) {
        for (int e = s; e <= maxSmall; ++e) {
            dp[e] = addmod(dp[e], dp[e - s]);
        }
    }

    long long res = 1;
    for (const auto& pr : small) {
        int e = pr.first;
        int cnt = pr.second;
        int val = dp[e];
        if (val == 0) return 0;
        res = (res * mod_pow(val, cnt)) % MOD;
    }

    if (big.empty()) {
        return static_cast<int>(res * entry.weight_mod % MOD);
    }

    array<int, K + 1> C{};
    C[0] = 1;
    int deg = 0;
    for (int s : entry.sums) {
        for (int i = deg; i >= 0; --i) {
            C[i + s] = submod(C[i + s], C[i]);
        }
        deg += s;
    }

    Vec rec{};
    for (int i = 1; i <= K; ++i) {
        rec[i - 1] = submod(0, C[i]);
    }

    const int* init = dp.data();
    for (const auto& pr : big) {
        int e = pr.first;
        int cnt = pr.second;
        int val = linear_rec(rec, init, e);
        if (val == 0) return 0;
        res = (res * mod_pow(val, cnt)) % MOD;
    }

    return static_cast<int>(res * entry.weight_mod % MOD);
}

static int compute_F(int n, unsigned threads) {
    ExpCounts counts = build_exp_counts(n);
    const auto& entries = get_entries();

    vector<int> active;
    active.reserve(entries.size());
    if (counts.hasExp1) {
        for (size_t i = 0; i < entries.size(); ++i) {
            if (entries[i].gcd == 1) active.push_back(static_cast<int>(i));
        }
    } else {
        for (size_t i = 0; i < entries.size(); ++i) {
            active.push_back(static_cast<int>(i));
        }
    }

    if (active.empty()) return 0;

    unsigned tcount = threads == 0 ? 1u : threads;
    if (tcount > active.size()) tcount = static_cast<unsigned>(active.size());

    const int inv288 = mod_pow(288, MOD - 2);

    if (tcount <= 1 || active.size() < 64) {
        vector<int> dp(counts.maxSmall + 1);
        long long total = 0;
        for (int idx : active) {
            total += compute_contribution(entries[idx], counts.small, counts.big, counts.maxSmall, dp);
            if (total >= MOD) total -= MOD;
        }
        return static_cast<int>(total * inv288 % MOD);
    }

    atomic<size_t> next(0);
    vector<long long> partial(tcount, 0);
    vector<thread> pool;
    pool.reserve(tcount);

    for (unsigned t = 0; t < tcount; ++t) {
        pool.emplace_back([&, t]() {
            vector<int> dp(counts.maxSmall + 1);
            long long local = 0;
            while (true) {
                size_t pos = next.fetch_add(1);
                if (pos >= active.size()) break;
                int idx = active[pos];
                local += compute_contribution(entries[idx], counts.small, counts.big, counts.maxSmall, dp);
                if (local >= MOD) local -= MOD;
            }
            partial[t] = local % MOD;
        });
    }

    for (auto& th : pool) th.join();

    long long total = 0;
    for (long long v : partial) {
        total += v;
        if (total >= MOD) total -= MOD;
    }

    return static_cast<int>(total * inv288 % MOD);
}

static bool run_validation(unsigned threads) {
    struct Case {
        int n;
        int expected;
    };

    const Case cases[] = {
        {25, 4933},
        {100, 693952493},
        {1000, 6364496},
    };

    bool ok = true;
    for (const auto& c : cases) {
        int got = compute_F(c.n, threads);
        if (got != c.expected) {
            cerr << "Validation failed for n=" << c.n
                 << ": got " << got << ", expected " << c.expected << "\n";
            ok = false;
        }
    }

    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 1'000'000;
    unsigned threads = thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = min(threads, 8u);
    bool validate = true;

    // Optional CLI: ./a.out [n] [threads] [validate(0/1)]
    if (argc >= 2) n = stoi(argv[1]);
    if (argc >= 3) threads = static_cast<unsigned>(stoul(argv[2]));
    if (argc >= 4) validate = (stoi(argv[3]) != 0);

    if (threads == 0) threads = 1;

    if (validate && !run_validation(min(threads, 4u))) {
        return 1;
    }

    cout << compute_F(n, threads) << "\n";
    return 0;
}
