#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;

using u64 = unsigned long long;
using u128 = unsigned __int128;

static u64 gcd_u64(u64 a, u64 b) {
    while (b) {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static u64 mod_inv(u64 a, u64 mod) {
    // mod assumed > 1 and gcd(a,mod)=1.
    __int128 t = 0, newt = 1;
    __int128 r = (u64)mod, newr = (u64)a;
    while (newr != 0) {
        u64 q = (u64)(r / newr);
        __int128 tmp_t = t - (__int128)q * newt;
        t = newt;
        newt = tmp_t;
        __int128 tmp_r = r - (__int128)q * newr;
        r = newr;
        newr = tmp_r;
    }
    if (r != 1) return 0;
    t %= (__int128)mod;
    if (t < 0) t += (__int128)mod;
    return (u64)t;
}

struct CRTRes {
    u64 r;  // residue in [0, m)
    u64 m;  // modulus
    bool ok;
};

static CRTRes crt_merge(u64 r1, u64 m1, u64 r2, u64 m2, u64 maxN) {
    // Solve x ≡ r1 (mod m1), x ≡ r2 (mod m2)
    if (m1 == 1) return {r2 % m2, m2, true};
    if (m2 == 1) return {r1 % m1, m1, true};
    u64 g = gcd_u64(m1, m2);
    if ((r1 % g) != (r2 % g)) return {0, 0, false};

    u64 m1_g = m1 / g;
    u64 m2_g = m2 / g;

    // Compute t = ((r2-r1)/g) * inv(m1/g) mod (m2/g)
    long long diff = (long long)r2 - (long long)r1; // divisible by g
    long long diff_g = diff / (long long)g;

    u64 inv = mod_inv(m1_g % m2_g, m2_g);
    u64 t = (u64)(((__int128)diff_g % (long long)m2_g + (long long)m2_g) % (long long)m2_g);
    t = (u64)((u128)t * inv % m2_g);

    u128 newm128 = (u128)m1_g * (u128)m2; // m1/g * m2
    u128 newr128 = (u128)r1 + (u128)m1 * (u128)t;

    // If modulus exceeds the numeric search range, only the smallest nonnegative solution can matter.
    const u64 clamp_mod = maxN + 1;
    if (newm128 > (u128)clamp_mod) {
        u128 rr = newr128 % newm128;
        if (rr > (u128)maxN) return {0, 0, false};
        return {(u64)rr, clamp_mod, true};
    }
    u64 newm = (u64)newm128;
    u64 newr = (u64)(newr128 % newm128);
    return {newr, newm, true};
}

static vector<int> sieve_primes(int n) {
    vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; (long long)i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    vector<int> primes;
    for (int i = 2; i <= n; ++i) if (is_prime[i]) primes.push_back(i);
    return primes;
}

static bool factor_squarefree_bounded(u64 x, int bound_exclusive, const vector<int> &primes,
                                     vector<int> &factors_out) {
    // Factor x by trial division using primes list.
    // Return true iff x is squarefree, all prime factors are < bound_exclusive, and x is odd.
    factors_out.clear();
    if (x == 0) return false;
    if ((x & 1ULL) == 0) return false; // even not allowed (except overall m=2, handled separately)

    u64 n = x;
    for (int p : primes) {
        if ((u64)p * (u64)p > n) break;
        if (p >= bound_exclusive) break;
        if (n % (u64)p == 0) {
            factors_out.push_back(p);
            n /= (u64)p;
            if (n % (u64)p == 0) return false; // not squarefree
            while (n % (u64)p == 0) n /= (u64)p; // shouldn't happen, but keep safe
        }
    }
    if (n > 1) {
        if (n >= (u64)bound_exclusive) return false;
        // n must be prime at this point (or a product of >= bound primes, which would be caught by n>=bound if bound small)
        factors_out.push_back((int)n);
    }
    return true;
}

static bool check_solution(u64 m, const vector<int> &prime_factors) {
    // prime_factors are distinct and should be the full factorization of m.
    u64 mp3 = m + 3;
    for (int p : prime_factors) {
        int d = p - 1;
        if (d == 0) return false;
        if (mp3 % (u64)d != 0) return false;
    }
    return true;
}

struct Solver {
    u64 N;
    vector<int> primes;      // primes up to sqrt(N)+4
    vector<int> primes_small; // primes up to 1e6 (same list reused)
    unordered_set<u64> found;

    static constexpr u64 MAX_ENUM = 10000; // target max candidates per branch

    Solver(u64 N_) : N(N_) {}

    pair<u64,u64> congruence_for_remaining(const vector<int> &B, u64 P) {
        // Compute R ≡ r (mod M) implied by primes in B for remaining product R.
        u64 r = 0;
        u64 M = 1;
        for (int p : B) {
            u64 mod = (u64)(p - 1);
            if (mod == 0) return {0, 0};
            u64 a = (P / (u64)p) % mod;
            u64 g = gcd_u64(a, mod);
            if (3 % g != 0) return {0, 0};
            u64 mod2 = mod / g;
            if (mod2 == 1) {
                // no constraint
                continue;
            }
            u64 a2 = a / g;
            u64 rhs = (mod2 - ((u64)(3 / g) % mod2)) % mod2; // -3/g mod mod2
            u64 inv = mod_inv(a2 % mod2, mod2);
            if (inv == 0) return {0, 0};
            u64 rp = (u64)((u128)rhs * inv % mod2);

            CRTRes merged = crt_merge(r, M, rp, mod2, N);
            if (!merged.ok) return {0, 0};
            r = merged.r;
            M = merged.m;
        }
        return {r, M};
    }

    u64 max_remaining_product(int bound_exclusive, u64 cap) const {
        // Maximum squarefree odd product using primes < bound_exclusive, capped at cap.
        u128 prod = 1;
        for (int p : primes) {
            if (p == 2) continue;
            if (p >= bound_exclusive) break;
            if (prod > (u128)cap / (u64)p) return cap;
            prod *= (u64)p;
        }
        if (prod > (u128)cap) return cap;
        return (u64)prod;
    }

    void enumerate_and_check(const vector<int> &B, u64 P, int min_big, u64 r, u64 M) {
        u64 limit = N / P;
        limit = min(limit, max_remaining_product(min_big, limit));
        if (M == 0) return;
        // Enumerate R = r + t*M
        vector<int> factorsR;
        vector<int> factorsM;
        factorsM.reserve(B.size() + 16);

        for (u64 R = r; R <= limit; ) {
            if (R != 0) {
                if (factor_squarefree_bounded(R, min_big, primes_small, factorsR)) {
                    // Merge factors: B (largest primes) + factorsR
                    factorsM.clear();
                    for (int p : B) factorsM.push_back(p);
                    for (int p : factorsR) factorsM.push_back(p);

                    u64 m = P * R;
                    if (m <= N) {
                        // Ensure B primes are indeed >= min_big and remaining primes < min_big already checked.
                        if (check_solution(m, factorsM)) {
                            found.insert(m);
                        }
                    }
                }
            }
            // next
            if (R > limit - M) break;
            R += M;
        }
    }

    void dfs(int max_prime_idx, vector<int> &B, u64 P) {
        int min_big = B.back();
        // Compute congruence for remaining product R
        auto [r, M] = congruence_for_remaining(B, P);
        if (M == 0) return;

        u64 limit = N / P;
        limit = min(limit, max_remaining_product(min_big, limit));
        if (limit == 0) return;
        if (r > limit) {
            // There might be a larger representative r + t*M, but r is in [0,M). If r>limit and M>limit, no candidates.
            // If M <= limit, then r + t*M could still be <= limit. Handle by shifting to first >= r.
            // We'll compute using floor division.
        }

        // Quick upper bound on count: limit / M + 1
        u64 approx_cnt = limit / M + 1;
        if (approx_cnt <= MAX_ENUM) {
            enumerate_and_check(B, P, min_big, r, M);
            return;
        }

        // Need to include more large primes (next prime factor) to reduce enumeration.
        // Choose next prime q < min_big.
        for (int idx = max_prime_idx - 1; idx >= 0; --idx) {
            int q = primes[idx];
            if (q == 2) continue; // even factors not allowed in composite solutions
            u128 P2 = (u128)P * (u64)q;
            if (P2 > (u128)N) continue;
            B.push_back(q);
            dfs(idx, B, (u64)P2);
            B.pop_back();
        }
    }

    u64 solve() {
        // Handle trivial solutions.
        found.clear();
        found.insert(1);
        if (N >= 2) found.insert(2);
        if (N >= 3) found.insert(3);
        if (N >= 5) found.insert(5);

        int prime_limit = (int)(sqrt((long double)N) + 4);
        primes = sieve_primes(prime_limit);
        primes_small = primes; // same for trial division

        // Exclude 2 from being a largest prime factor in composite enumeration.
        // Iterate over odd primes >= 7 (since we already inserted 3 and 5).
        for (int i = (int)primes.size() - 1; i >= 0; --i) {
            int pmax = primes[i];
            if (pmax < 7) break;
            vector<int> B;
            B.push_back(pmax);
            dfs(i, B, (u64)pmax);
        }

        u128 sum = 0;
        for (u64 m : found) {
            if (m <= N) sum += (u128)m;
        }
        return (u64)sum;
    }
};

static u64 brute_sum(u64 N) {
    // Brute using SPF up to N (N should be small, <= 1e7 here).
    int n = (int)N;
    vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) spf[i] = i;
    for (int i = 2; (long long)i * i <= n; ++i) {
        if (spf[i] != i) continue;
        for (long long j = 1LL * i * i; j <= n; j += i) {
            if (spf[(int)j] == (int)j) spf[(int)j] = i;
        }
    }

    u128 sum = 0;
    for (int m = 1; m <= n; ++m) {
        int x = m;
        bool ok = true;
        while (x > 1) {
            int p = spf[x];
            int cnt = 0;
            while (x % p == 0) {
                x /= p;
                if (++cnt > 1) {
                    ok = false;
                    break;
                }
            }
            if (!ok) break;
            if ((u64)(m + 3) % (u64)(p - 1) != 0) {
                ok = false;
                break;
            }
        }
        if (ok) sum += (u128)m;
    }
    return (u64)sum;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Validation checkpoints from the problem statement.
    {
        u64 s100 = brute_sum(100);
        if (s100 != 32) {
            cerr << "Checkpoint failed: S(100)=" << s100 << " (expected 32)\n";
            return 1;
        }
        u64 s1e6 = brute_sum(1000000ULL);
        if (s1e6 != 22868117ULL) {
            cerr << "Checkpoint failed: S(10^6)=" << s1e6 << " (expected 22868117)\n";
            return 1;
        }
    }

    const u64 N = 1000000000000ULL;
    Solver solver(N);
    u64 ans = solver.solve();
    cout << ans << "\n";
    return 0;
}
