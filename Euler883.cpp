#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

using namespace std;

namespace {

// floor(sqrt(x)) for 0 <= x <= 2^64-1.
uint64_t isqrt_u64(uint64_t x) {
    long double r = sqrtl(static_cast<long double>(x));
    uint64_t y = static_cast<uint64_t>(r);
    while (static_cast<__int128>(y + 1) * static_cast<__int128>(y + 1) <= x) {
        ++y;
    }
    while (static_cast<__int128>(y) * static_cast<__int128>(y) > x) {
        --y;
    }
    return y;
}

int64_t floor_div2(int64_t a) {
    if (a >= 0) return a >> 1;
    return -(((-a) + 1) >> 1);
}

int64_t ceil_div2(int64_t a) {
    if (a >= 0) return (a + 1) >> 1;
    return -(((-a)) >> 1);
}

int64_t mod_pos(int64_t a, int64_t m) {
    int64_t r = a % m;
    if (r < 0) r += m;
    return r;
}

// Count lattice points (u, v) with u^2 + u v + v^2 <= M.
// Returns (A, Aeq), where A counts all points (including origin),
// and Aeq counts points with u ≡ v (mod 3) (including origin).
pair<uint64_t, uint64_t> count_both(uint64_t M) {
    const uint64_t fourM = 4ULL * M;
    const uint64_t umax = isqrt_u64(fourM / 3ULL);

    uint64_t all = 0;
    uint64_t eq = 0;

    for (uint64_t u = 0; u <= umax; ++u) {
        const uint64_t D = fourM - 3ULL * u * u;
        const uint64_t s = isqrt_u64(D);

        const int64_t uu = static_cast<int64_t>(u);
        const int64_t ss = static_cast<int64_t>(s);
        const int64_t vmin = ceil_div2(-uu - ss);
        const int64_t vmax = floor_div2(-uu + ss);

        uint64_t len = 0;
        if (vmax >= vmin) len = static_cast<uint64_t>(vmax - vmin + 1);

        uint64_t cong = 0;
        if (len) {
            const int r = static_cast<int>(u % 3ULL);
            const int64_t first = vmin + mod_pos(static_cast<int64_t>(r) - vmin, 3);
            if (first <= vmax) {
                cong = static_cast<uint64_t>((vmax - first) / 3 + 1);
            }
        }

        if (u == 0) {
            all += len;
            eq += cong;
        } else {
            all += 2ULL * len;
            eq += 2ULL * cong;
        }
    }
    return {all, eq};
}

// Linear sieve for smallest prime factor and Mobius function.
void linear_sieve_spf_mu(int nmax, vector<int>& spf, vector<int8_t>& mu) {
    spf.assign(nmax + 1, 0);
    mu.assign(nmax + 1, 0);
    vector<int> primes;
    primes.reserve(nmax / 10);

    spf[1] = 1;
    mu[1] = 1;

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

// Compute f(n) used in the convolution.
int32_t compute_f_single(int n, const vector<int>& spf) {
    int x = n;
    int e3 = 0;
    while (x % 3 == 0) {
        x /= 3;
        ++e3;
    }

    long long tau_h2 = 1;
    long long prod1 = 1;
    int chi = 1;

    while (x > 1) {
        const int p = spf[x];
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        const long long factor = 2LL * e + 1;
        tau_h2 *= factor;

        const int pm3 = p % 3;
        if (pm3 == 1) {
            prod1 *= factor;
        } else if (pm3 == 2) {
            if (e & 1) chi = -chi;
        }
    }

    long long res;
    if (e3 > 0) {
        res = tau_h2 * (2LL * e3 - 1);
    } else {
        res = (tau_h2 + 1LL * chi * prod1) / 2LL;
    }
    return static_cast<int32_t>(res);
}

struct Result {
    int64_t T = 0;
    int64_t S = 0;
    uint64_t points = 0;
};

Result solve_R2(int64_t R2, unsigned numThreads, bool doCheckpoints) {
    const uint64_t N = static_cast<uint64_t>(R2) * static_cast<uint64_t>(R2);
    const int nmax = static_cast<int>(3 * R2);

    vector<int> spf;
    vector<int8_t> mu;
    linear_sieve_spf_mu(nmax, spf, mu);

    vector<int32_t> f(nmax + 1, 0);
    for (int n = 1; n <= nmax; ++n) {
        f[n] = compute_f_single(n, spf);
    }

    // h = f * mu (Dirichlet convolution).
    vector<int32_t> h(nmax + 1, 0);
    for (int k = 1; k <= nmax; ++k) {
        const int8_t muk = mu[k];
        if (muk == 0) continue;
        for (int d = 1, m = k; m <= nmax; ++d, m += k) {
            h[m] += static_cast<int32_t>(muk * f[d]);
        }
    }

    vector<uint64_t> F_base(R2 + 1, 0);
    vector<uint64_t> F_eq(R2 + 1, 0);

    atomic<int> nextT(1);
    const int chunk = 64;

    auto worker = [&]() {
        while (true) {
            const int start = nextT.fetch_add(chunk);
            if (start > R2) break;
            const int end = min<int>(R2, start + chunk - 1);
            for (int t = start; t <= end; ++t) {
                const uint64_t denom = static_cast<uint64_t>(t) * static_cast<uint64_t>(t);
                const uint64_t M = N / denom;
                auto [all, eq] = count_both(M);
                F_base[t] = all - 1ULL;
                F_eq[t] = eq - 1ULL;
            }
        }
    };

    if (numThreads == 0) numThreads = 1;
    vector<thread> threads;
    threads.reserve(numThreads);
    for (unsigned i = 0; i < numThreads; ++i) {
        threads.emplace_back(worker);
    }
    for (auto& th : threads) {
        th.join();
    }

    const uint64_t totalPoints = F_base[1];

    if (doCheckpoints) {
        if (totalPoints % 3ULL != 0ULL) {
            throw runtime_error("Validation failed: totalPoints not divisible by 3");
        }
        for (int t : {1, 2, 3, 10, static_cast<int>(R2)}) {
            if (t <= 0 || t > R2) continue;
            if (F_eq[t] > F_base[t]) {
                throw runtime_error("Validation failed: F_eq[t] > F_base[t]");
            }
        }
    }

    const uint64_t E = totalPoints / 3ULL;

    __int128 S = 0;
    for (int n = 1; n <= nmax; ++n) {
        const int32_t hn = h[n];
        if (hn == 0) continue;

        uint64_t Fn = 0;
        if (n % 3 != 0) {
            if (n <= R2) Fn = F_base[n];
        } else {
            const int t = n / 3;
            Fn = F_eq[t];
        }
        S += static_cast<__int128>(hn) * static_cast<__int128>(Fn);
    }

    const __int128 T = S - 2 * static_cast<__int128>(E);

    Result res;
    res.T = static_cast<int64_t>(T);
    res.S = static_cast<int64_t>(S);
    res.points = totalPoints;
    return res;
}

void run_small_validations(unsigned threads) {
    struct Case { int R2; int64_t expected; };
    const vector<Case> cases = {
        {1, 2},
        {4, 44},
        {20, 1302},
    };

    for (const auto& c : cases) {
        const Result r = solve_R2(c.R2, threads, true);
        if (r.T != c.expected) {
            cerr << "Validation failed for R2=" << c.R2
                 << " (r=" << (c.R2 / 2.0) << "): got " << r.T
                 << ", expected " << c.expected << "\n";
            exit(1);
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int64_t R2 = 2'000'000;
    unsigned threads = thread::hardware_concurrency();
    bool validate = true;

    // Optional CLI: ./a.out [R2] [threads] [validate(0/1)]
    if (argc >= 2) R2 = stoll(argv[1]);
    if (argc >= 3) threads = static_cast<unsigned>(stoul(argv[2]));
    if (argc >= 4) validate = (stoi(argv[3]) != 0);

    if (validate) {
        run_small_validations(min(threads, 4u));
    }

    const Result ans = solve_R2(R2, threads, true);
    cout << ans.T << "\n";
    return 0;
}
