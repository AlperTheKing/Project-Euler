#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr int64_t MOD = 214358881LL;     // 11^8
constexpr int64_t PHI = 194871710LL;     // 10 * 11^7
constexpr int64_t PHI2 = 389743420LL;    // 2 * PHI

struct Pow2Mod {
    int64_t mod;
    array<int64_t, 32> table{};
    explicit Pow2Mod(int64_t mod_) : mod(mod_) {
        table[0] = 2 % mod;
        for (size_t i = 1; i < table.size(); ++i) {
            table[i] = (table[i - 1] * table[i - 1]) % mod;
        }
    }

    int64_t pow(int64_t exp) const {
        int64_t res = 1;
        size_t bit = 0;
        while (exp > 0) {
            if (exp & 1LL) res = (res * table[bit]) % mod;
            exp >>= 1;
            ++bit;
        }
        return res;
    }
};

void sieve_mu_phi(int n, vector<int>& mu, vector<int>& phi) {
    mu.assign(n + 1, 0);
    phi.assign(n + 1, 0);
    vector<int> primes;
    vector<uint8_t> is_comp(n + 1, 0);

    mu[1] = 1;
    phi[1] = 1;

    for (int i = 2; i <= n; ++i) {
        if (!is_comp[i]) {
            primes.push_back(i);
            mu[i] = -1;
            phi[i] = i - 1;
        }
        for (int p : primes) {
            int64_t v = 1LL * i * p;
            if (v > n) break;
            is_comp[v] = 1;
            if (i % p == 0) {
                mu[v] = 0;
                phi[v] = phi[i] * p;
                break;
            }
            mu[v] = -mu[i];
            phi[v] = phi[i] * (p - 1);
        }
    }
}

int64_t count_M_mod2phi(int n, const vector<int>& mu, const vector<int>& n_div) {
    int64_t total = 0;
    int64_t mertens = 0;
    for (int d = 1; d <= n; ++d) {
        int mu_d = mu[d];
        if (mu_d == 0) continue;
        int64_t k = n_div[d];
        __int128 cube = (k + 1);
        cube = cube * cube * cube;
        int64_t term = static_cast<int64_t>(cube % PHI2);
        total += static_cast<int64_t>(mu_d) * term;
        mertens += mu_d;
        if (total >= PHI2 || total <= -PHI2) total %= PHI2;
        if (mertens >= PHI2 || mertens <= -PHI2) mertens %= PHI2;
    }
    total = (total - mertens) % PHI2;
    if (total < 0) total += PHI2;
    return total;
}

int64_t compute_E(int n, unsigned threads) {
    vector<int> mu;
    vector<int> phi;
    sieve_mu_phi(n, mu, phi);

    vector<int> n_div(n + 1, 0);
    for (int d = 1; d <= n; ++d) n_div[d] = n / d;

    int64_t sum_phi_mod = 0;
    for (int i = 1; i <= n; ++i) {
        sum_phi_mod += phi[i];
        if (sum_phi_mod >= MOD) sum_phi_mod -= MOD;
    }
    int64_t D_mod = (6LL * sum_phi_mod) % MOD;

    int64_t N_mod2phi = count_M_mod2phi(n, mu, n_div);
    int64_t N_mod_phi = N_mod2phi % PHI;
    int64_t Nprime_mod2phi = (N_mod2phi + PHI2 - 1) % PHI2;
    if (Nprime_mod2phi & 1LL) {
        cerr << "Validation failed: N' not even." << '\n';
        return -1;
    }
    int64_t Nhalf_mod_phi = (Nprime_mod2phi / 2) % PHI;

    Pow2Mod pow2(MOD);
    int64_t pow2_N = pow2.pow(N_mod_phi);
    int64_t pow2_half = pow2.pow(Nhalf_mod_phi);

    if (threads == 0) threads = 1;
    if (n < 100000) threads = 1;
    threads = min<unsigned>(threads, static_cast<unsigned>(n));

    vector<int64_t> partial(threads, 0);
    atomic<int> nextL{1};
    const int chunk = 64;

    auto worker = [&](unsigned idx) {
        int64_t local = 0;
        while (true) {
            int start = nextL.fetch_add(chunk, memory_order_relaxed);
            if (start > n) break;
            int end = min(n, start + chunk - 1);

            for (int L = start; L <= end; ++L) {
                int m = n / L;
                int64_t total = 1;
                for (int d = 1; d <= m; ++d) {
                    int mu_d = mu[d];
                    if (mu_d == 0) continue;
                    int md = m / d;
                    int64_t term = 1LL * md * n_div[d] - 1LL * L * md * (md + 1) / 2;
                    total += static_cast<int64_t>(mu_d) * term;
                }
                int64_t g_mod = total % PHI;
                if (g_mod < 0) g_mod += PHI;
                int64_t exp = Nhalf_mod_phi - g_mod;
                if (exp < 0) exp += PHI;
                int64_t pow = pow2.pow(exp);
                int64_t mL = (6LL * phi[L]) % MOD;
                local += (mL * pow) % MOD;
                if (local >= MOD) local %= MOD;
            }
        }
        partial[idx] = local % MOD;
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    int64_t S = 0;
    for (int64_t v : partial) {
        S += v;
        if (S >= MOD) S %= MOD;
    }

    int64_t F = (1 + (pow2_half * D_mod) % MOD - S) % MOD;
    if (F < 0) F += MOD;
    int64_t E = (pow2_N - F) % MOD;
    if (E < 0) E += MOD;
    return E;
}

bool validate() {
    struct Test {
        int n;
        int64_t expected;
    };
    const Test tests[] = {
        {1, 103},
        {2, 520447},
        {10, 82608406},
        {500, 13801403},
    };

    for (const auto& test : tests) {
        int64_t got = compute_E(test.n, 1);
        if (got != test.expected) {
            cerr << "Validation failed for n=" << test.n
                 << ": got " << got << ", expected " << test.expected << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (!validate()) return 1;

    int n = 10000000;
    if (argc > 1) {
        n = max(1, atoi(argv[1]));
    }
    unsigned threads = thread::hardware_concurrency();
    int64_t result = compute_E(n, threads);
    if (result < 0) return 1;

    cout << result << '\n';
    return 0;
}
