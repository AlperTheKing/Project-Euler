#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr uint64_t kTargetN = 1'000'000'000'000'000'000ULL;
constexpr uint64_t kMod = 1'000'000'000ULL;

int floor_log2_u64(uint64_t x) {
    int r = 0;
    while (x > 1) {
        x >>= 1;
        ++r;
    }
    return r;
}

bool pow_leq(uint64_t base, int exp, uint64_t limit) {
    unsigned __int128 res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= base;
        if (res > limit) return false;
    }
    return true;
}

uint64_t int_nth_root(uint64_t n, int k) {
    if (k == 1 || n <= 1) return n;
    long double approx = powl(static_cast<long double>(n), 1.0L / k);
    uint64_t r = static_cast<uint64_t>(approx);
    if (r < 1) r = 1;
    while (pow_leq(r + 1, k, n)) ++r;
    while (!pow_leq(r, k, n)) --r;
    return r;
}

void mobius_sieve(int nmax, vector<int8_t>& mu) {
    mu.assign(nmax + 1, 0);
    vector<int> primes;
    primes.reserve(nmax / 2);
    mu[1] = 1;
    vector<int> spf(nmax + 1, 0);
    for (int i = 2; i <= nmax; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            long long v = 1LL * p * i;
            if (v > nmax) break;
            spf[static_cast<int>(v)] = p;
            if (i % p == 0) {
                mu[static_cast<int>(v)] = 0;
                break;
            }
            mu[static_cast<int>(v)] = static_cast<int8_t>(-mu[i]);
        }
    }
}

// Count integers r in [2, M] that are not perfect powers.
uint64_t power_free_count_upto(uint64_t M, const vector<int8_t>& mu) {
    if (M < 2) return 0;
    int L = floor_log2_u64(M);
    int max_d = min(L, static_cast<int>(mu.size()) - 1);
    int64_t total = 0;
    for (int d = 1; d <= max_d; ++d) {
        const int8_t md = mu[d];
        if (md == 0) continue;
        uint64_t root = int_nth_root(M, d);
        if (root >= 2) {
            total += static_cast<int64_t>(md) * static_cast<int64_t>(root - 1);
        }
    }
    return static_cast<uint64_t>(total);
}

vector<vector<int64_t>> build_coprime_table(int L, const vector<int8_t>& mu) {
    vector<vector<int64_t>> F(L + 1, vector<int64_t>(L + 1, 0));
    for (int A = 1; A <= L; ++A) {
        for (int B = 1; B <= L; ++B) {
            int m = min(A, B);
            int64_t total = 0;
            for (int d = 1; d <= m; ++d) {
                int8_t md = mu[d];
                if (md == 0) continue;
                total += static_cast<int64_t>(md) * (A / d) * (B / d);
            }
            F[A][B] = total;
        }
    }
    return F;
}

vector<uint64_t> build_countA(uint64_t N, const vector<int8_t>& mu, int L) {
    vector<uint64_t> roots(L + 2, 0);
    for (int k = 1; k <= L; ++k) {
        roots[k] = int_nth_root(N, k);
    }
    roots[L + 1] = 1;

    vector<uint64_t> pf(L + 2, 0);
    for (int k = 1; k <= L + 1; ++k) {
        pf[k] = power_free_count_upto(roots[k], mu);
    }

    vector<uint64_t> countA(L + 1, 0);
    for (int A = 1; A <= L; ++A) {
        countA[A] = pf[A] - pf[A + 1];
    }
    return countA;
}

uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t mod) {
    return static_cast<uint64_t>((static_cast<unsigned __int128>(a) * b) % mod);
}

unsigned __int128 compute_irrational_exact(const vector<uint64_t>& countA,
                                           const vector<vector<int64_t>>& coprime,
                                           unsigned threads) {
    const int L = static_cast<int>(countA.size()) - 1;
    if (threads <= 1 || L <= 1) {
        unsigned __int128 irr = 0;
        for (int A = 1; A <= L; ++A) {
            uint64_t ca = countA[A];
            if (ca == 0) continue;
            for (int B = 1; B <= L; ++B) {
                uint64_t cb = countA[B];
                if (cb == 0) continue;
                unsigned __int128 pairs = 0;
                if (A == B) {
                    if (ca < 2) continue;
                    pairs = static_cast<unsigned __int128>(ca) * (ca - 1);
                } else {
                    pairs = static_cast<unsigned __int128>(ca) * cb;
                }
                irr += pairs * static_cast<unsigned __int128>(coprime[A][B]);
            }
        }
        return irr;
    }

    threads = min<unsigned>(threads, static_cast<unsigned>(L));
    vector<unsigned __int128> partial(threads, 0);
    atomic<int> nextA(1);
    const int chunk = 2;

    auto worker = [&](unsigned idx) {
        unsigned __int128 local = 0;
        while (true) {
            int start = nextA.fetch_add(chunk);
            if (start > L) break;
            int end = min(L, start + chunk - 1);
            for (int A = start; A <= end; ++A) {
                uint64_t ca = countA[A];
                if (ca == 0) continue;
                for (int B = 1; B <= L; ++B) {
                    uint64_t cb = countA[B];
                    if (cb == 0) continue;
                    unsigned __int128 pairs = 0;
                    if (A == B) {
                        if (ca < 2) continue;
                        pairs = static_cast<unsigned __int128>(ca) * (ca - 1);
                    } else {
                        pairs = static_cast<unsigned __int128>(ca) * cb;
                    }
                    local += pairs * static_cast<unsigned __int128>(coprime[A][B]);
                }
            }
        }
        partial[idx] = local;
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for (auto& th : pool) th.join();

    unsigned __int128 total = 0;
    for (auto v : partial) total += v;
    return total;
}

uint64_t compute_irrational_mod(const vector<uint64_t>& countA,
                                const vector<vector<int64_t>>& coprime,
                                unsigned threads,
                                uint64_t mod) {
    const int L = static_cast<int>(countA.size()) - 1;
    if (threads == 0) threads = 1;
    if (threads <= 1 || L <= 1) {
        uint64_t irr = 0;
        for (int A = 1; A <= L; ++A) {
            uint64_t ca = countA[A];
            if (ca == 0) continue;
            uint64_t ca_mod = ca % mod;
            for (int B = 1; B <= L; ++B) {
                uint64_t cb = countA[B];
                if (cb == 0) continue;
                uint64_t pairs = 0;
                if (A == B) {
                    if (ca < 2) continue;
                    pairs = mod_mul(ca_mod, (ca - 1) % mod, mod);
                } else {
                    pairs = mod_mul(ca_mod, cb % mod, mod);
                }
                uint64_t term = mod_mul(pairs,
                                        static_cast<uint64_t>(coprime[A][B]) % mod,
                                        mod);
                irr += term;
                if (irr >= mod) irr %= mod;
            }
        }
        return irr % mod;
    }

    threads = min<unsigned>(threads, static_cast<unsigned>(L));
    vector<uint64_t> partial(threads, 0);
    atomic<int> nextA(1);
    const int chunk = 2;

    auto worker = [&](unsigned idx) {
        uint64_t local = 0;
        while (true) {
            int start = nextA.fetch_add(chunk);
            if (start > L) break;
            int end = min(L, start + chunk - 1);
            for (int A = start; A <= end; ++A) {
                uint64_t ca = countA[A];
                if (ca == 0) continue;
                uint64_t ca_mod = ca % mod;
                for (int B = 1; B <= L; ++B) {
                    uint64_t cb = countA[B];
                    if (cb == 0) continue;
                    uint64_t pairs = 0;
                    if (A == B) {
                        if (ca < 2) continue;
                        pairs = mod_mul(ca_mod, (ca - 1) % mod, mod);
                    } else {
                        pairs = mod_mul(ca_mod, cb % mod, mod);
                    }
                    uint64_t term = mod_mul(pairs,
                                            static_cast<uint64_t>(coprime[A][B]) % mod,
                                            mod);
                    local += term;
                    if (local >= mod) local %= mod;
                }
            }
        }
        partial[idx] = local % mod;
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for (auto& th : pool) th.join();

    uint64_t total = 0;
    for (uint64_t v : partial) {
        total += v;
        if (total >= mod) total %= mod;
    }
    return total % mod;
}

unsigned __int128 compute_D_exact(uint64_t N, unsigned threads) {
    const int L = floor_log2_u64(N);
    vector<int8_t> mu;
    mobius_sieve(L, mu);
    auto coprime = build_coprime_table(L, mu);
    auto countA = build_countA(N, mu, L);

    unsigned __int128 irrational = compute_irrational_exact(countA, coprime, threads);
    unsigned __int128 rational = static_cast<unsigned __int128>(coprime[L][L]);
    return rational + irrational;
}

uint64_t compute_D_mod(uint64_t N, unsigned threads, uint64_t mod) {
    const int L = floor_log2_u64(N);
    vector<int8_t> mu;
    mobius_sieve(L, mu);
    auto coprime = build_coprime_table(L, mu);
    auto countA = build_countA(N, mu, L);

    uint64_t irrational = compute_irrational_mod(countA, coprime, threads, mod);
    uint64_t rational = static_cast<uint64_t>(coprime[L][L]) % mod;
    return (rational + irrational) % mod;
}

string to_string_u128(unsigned __int128 x) {
    if (x == 0) return "0";
    string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}

bool run_validations(unsigned threads) {
    struct TestCase {
        uint64_t n;
        uint64_t expected;
    };
    const TestCase tests[] = {
        {5, 13},
        {10, 69},
        {100, 9607},
        {10000, 99959605},
    };

    for (const auto& tc : tests) {
        unsigned __int128 value = compute_D_exact(tc.n, threads);
        uint64_t got = static_cast<uint64_t>(value);
        if (got != tc.expected) {
            cerr << "Validation failed: D(" << tc.n << ") = " << got
                 << " (expected " << tc.expected << ")\n";
            return false;
        }
    }
    cerr << "Validation checkpoints passed.\n";
    return true;
}

} // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    uint64_t N = kTargetN;
    unsigned threads = thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    bool validate = true;

    if (argc >= 2) N = stoull(argv[1]);
    if (argc >= 3) threads = max(1u, static_cast<unsigned>(stoul(argv[2])));
    if (argc >= 4) validate = (stoi(argv[3]) != 0);

    if (validate && !run_validations(min(threads, 4u))) {
        return 1;
    }

    if (N == kTargetN) {
        uint64_t answer = compute_D_mod(N, threads, kMod);
        cout << setfill('0') << setw(9) << answer << "\n";
    } else {
        unsigned __int128 exact = compute_D_exact(N, threads);
        cout << to_string_u128(exact) << "\n";
    }
    return 0;
}
