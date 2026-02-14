#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <functional>

using namespace std;

namespace {

constexpr int64_t MOD = 1'000'000'007LL;

int64_t mod_norm(__int128 v) {
    v %= MOD;
    if (v < 0) v += MOD;
    return static_cast<int64_t>(v);
}

int64_t count_intervals(int terms) {
    if (terms < 5) return 0;
    __int128 a = terms - 4;
    __int128 b = terms - 3;
    return mod_norm(a * b / 2);
}

uint64_t pow4(uint64_t x) {
    __int128 v = static_cast<__int128>(x) * x;
    v *= v;
    return static_cast<uint64_t>(v);
}

uint64_t floor_root4(uint64_t n) {
    long double approx = pow(static_cast<long double>(n), 0.25L);
    uint64_t x = static_cast<uint64_t>(approx);
    while (pow4(x + 1) <= n) ++x;
    while (pow4(x) > n) --x;
    return x;
}

vector<int64_t> sieve_phi(int n) {
    vector<int64_t> phi(n + 1);
    for (int i = 0; i <= n; ++i) phi[i] = i;
    for (int i = 2; i <= n; ++i) {
        if (phi[i] == i) {
            for (int j = i; j <= n; j += i) {
                phi[j] -= phi[j] / i;
            }
        }
    }
    return phi;
}

int64_t count_geometric_range(uint64_t n, const vector<int64_t>& phi,
                              uint64_t start, uint64_t end) {
    int64_t sum = 0;
    for (uint64_t p = start; p <= end; ++p) {
        __int128 pk = static_cast<__int128>(p) * p;
        pk *= pk;  // p^4
        while (pk <= n) {
            uint64_t div = n / static_cast<uint64_t>(pk);
            __int128 term = static_cast<__int128>(phi[p]) * div;
            sum += mod_norm(term);
            if (sum >= MOD) sum -= MOD;
            pk *= p;
        }
    }
    return sum;
}

int64_t count_geometric(uint64_t n, uint64_t root4) {
    if (root4 < 2) return 0;
    vector<int64_t> phi = sieve_phi(static_cast<int>(root4));

    int thread_count = static_cast<int>(thread::hardware_concurrency());
    if (thread_count <= 0) thread_count = 1;
    if (root4 < 4000) thread_count = 1;
    thread_count = min(thread_count, static_cast<int>(root4 - 1));

    vector<int64_t> partial(thread_count, 0);
    vector<thread> workers;
    workers.reserve(thread_count);

    for (int t = 0; t < thread_count; ++t) {
        uint64_t start = 2 + (root4 - 1) * t / thread_count;
        uint64_t end = 2 + (root4 - 1) * (t + 1) / thread_count - 1;
        workers.emplace_back([&, t, start, end]() {
            partial[t] = count_geometric_range(n, phi, start, end);
        });
    }
    for (auto& th : workers) th.join();

    int64_t total = 0;
    for (int64_t v : partial) {
        total += v;
        if (total >= MOD) total -= MOD;
    }
    return total;
}

int count_terms(uint64_t n, uint64_t a0, uint64_t a1, uint64_t m, bool plus) {
    int terms = 2;
    for (;;) {
        __int128 a2 = plus ? static_cast<__int128>(m) * a1 + a0
                           : static_cast<__int128>(m) * a1 - a0;
        if (a2 <= a1 || a2 > n) break;
        ++terms;
        a0 = a1;
        a1 = static_cast<uint64_t>(a2);
    }
    return terms;
}

int64_t count_consecutive(uint64_t n) {
    if (n < 5) return 0;
    __int128 a = static_cast<__int128>(n - 4);
    __int128 b = static_cast<__int128>(n - 3);
    return mod_norm(a * b / 2);
}

int64_t count_hybrid(uint64_t n) {
    if (n < 54) return 0;
    int terms = 1;  // leading 1
    __int128 v = 2;
    while (v <= n) {
        ++terms;
        v *= 3;
    }
    if (terms < 5) return 0;
    return (terms - 4) % MOD;
}

int64_t count_sporadic(uint64_t n) {
    static const uint64_t last_terms[] = {
        6, 9, 9, 16, 12, 20, 48, 60, 9, 16
    };
    int64_t total = 0;
    for (uint64_t v : last_terms) {
        if (v <= n) ++total;
    }
    return total % MOD;
}

int64_t compute_G(uint64_t n) {
    if (n < 5) return 0;

    uint64_t root4 = floor_root4(n);
    int64_t ans = 0;

    ans = (ans + count_consecutive(n)) % MOD;
    ans = (ans + count_geometric(n, root4)) % MOD;

    // Minus recurrence: a_{k+1} = m a_k - a_{k-1}, k=1 sequences start (1,m).
    for (uint64_t m = 3; m <= root4; ++m) {
        int terms = count_terms(n, 1, m, m, false);
        ans = (ans + count_intervals(terms)) % MOD;
    }

    // Extra minus families.
    ans = (ans + count_intervals(count_terms(n, 1, 2, 3, false))) % MOD;
    ans = (ans + count_intervals(count_terms(n, 1, 3, 4, false))) % MOD;

    // Plus recurrence: a_{k+1} = m a_k + a_{k-1}, |k| = 1 families.
    ans = (ans + count_intervals(count_terms(n, 1, 2, 1, true))) % MOD;
    for (uint64_t m = 2; m <= root4; ++m) {
        int terms = count_terms(n, 1, m, m, true);
        ans = (ans + count_intervals(terms)) % MOD;
    }

    // Plus recurrence with |k| = 2 (m = 2, start 1,3).
    ans = (ans + count_intervals(count_terms(n, 1, 3, 2, true))) % MOD;

    ans = (ans + count_hybrid(n)) % MOD;
    ans = (ans + count_sporadic(n)) % MOD;

    return ans;
}

}  // namespace

int main() {
    const uint64_t N = 1'000'000'000'000'000'000ULL;

    const pair<uint64_t, int64_t> checks[] = {
        {6, 4},
        {10, 26},
        {100, 4710},
        {1000, 496805},
    };
    for (const auto& chk : checks) {
        int64_t got = compute_G(chk.first);
        if (got != chk.second) {
            cerr << "Validation failure: G(" << chk.first << ") = "
                 << got << ", expected " << chk.second << '\n';
            return 1;
        }
    }

    cout << compute_G(N) << '\n';
    return 0;
}
