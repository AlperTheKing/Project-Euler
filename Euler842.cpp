#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;

namespace {

constexpr int64_t MOD = 1'000'000'007LL;
using BigFloat = boost::multiprecision::cpp_dec_float_50;
const BigFloat SCALE_BF("1e18");

int64_t mod_mul(int64_t a, int64_t b) {
    return static_cast<int64_t>((__int128)a * b % MOD);
}

int64_t mod_pow(int64_t a, int64_t e) {
    int64_t r = 1 % MOD;
    int64_t x = (a % MOD + MOD) % MOD;
    while (e > 0) {
        if (e & 1) r = mod_mul(r, x);
        x = mod_mul(x, x);
        e >>= 1;
    }
    return r;
}

struct Precomp {
    vector<int64_t> fact;
    vector<int64_t> invfact;
    vector<int64_t> pow2;
    int64_t inv2;
};

Precomp precompute(int max_n) {
    Precomp p;
    p.fact.assign(max_n + 1, 1);
    p.invfact.assign(max_n + 1, 1);
    p.pow2.assign(max_n + 1, 1);
    for (int i = 1; i <= max_n; ++i) {
        p.fact[i] = mod_mul(p.fact[i - 1], i);
        p.pow2[i] = mod_mul(p.pow2[i - 1], 2);
    }
    p.invfact[max_n] = mod_pow(p.fact[max_n], MOD - 2);
    for (int i = max_n; i >= 1; --i) {
        p.invfact[i - 1] = mod_mul(p.invfact[i], i);
    }
    p.inv2 = mod_pow(2, MOD - 2);
    return p;
}

int64_t comb(int n, int k, const Precomp& p) {
    if (k < 0 || k > n) return 0;
    return mod_mul(p.fact[n], mod_mul(p.invfact[k], p.invfact[n - k]));
}

int64_t cycles_with_t_edges(int n, int t, const Precomp& p) {
    if (t == 0) return mod_mul(p.fact[n - 1], p.inv2);
    return mod_mul(p.fact[n - t - 1], p.pow2[t - 1]);
}

int64_t count_at_least_two(int n, int m, const Precomp& p) {
    int64_t total = cycles_with_t_edges(n, 0, p);
    int64_t none = 0;
    for (int t = 0; t <= m; ++t) {
        int64_t term = mod_mul(comb(m, t, p), cycles_with_t_edges(n, t, p));
        if (t & 1) {
            none = (none - term) % MOD;
        } else {
            none = (none + term) % MOD;
        }
    }
    if (none < 0) none += MOD;

    int64_t exactly_one = 0;
    for (int t = 0; t <= m - 1; ++t) {
        int64_t term = mod_mul(comb(m - 1, t, p), cycles_with_t_edges(n, t + 1, p));
        if (t & 1) {
            exactly_one = (exactly_one - term) % MOD;
        } else {
            exactly_one = (exactly_one + term) % MOD;
        }
    }
    if (exactly_one < 0) exactly_one += MOD;
    exactly_one = mod_mul(exactly_one, m % MOD);

    int64_t result = (total - none - exactly_one) % MOD;
    if (result < 0) result += MOD;
    return result;
}

struct PointKey {
    int64_t x;
    int64_t y;

    bool operator==(const PointKey& other) const noexcept {
        return x == other.x && y == other.y;
    }
};

struct PointKeyHash {
    size_t operator()(const PointKey& p) const noexcept {
        size_t h1 = std::hash<int64_t>{}(p.x);
        size_t h2 = std::hash<int64_t>{}(p.y);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

int64_t choose4(int n) {
    if (n < 4) return 0;
    return static_cast<int64_t>(n) * (n - 1) * (n - 2) * (n - 3) / 24;
}

int64_t round_to_int64(const BigFloat& v) {
    BigFloat scaled = v * SCALE_BF;
    if (scaled >= 0) {
        return static_cast<int64_t>(scaled + BigFloat("0.5"));
    }
    return -static_cast<int64_t>(-scaled + BigFloat("0.5"));
}

vector<int64_t> intersection_multiplicities(int n) {
    using boost::multiprecision::acos;
    using boost::multiprecision::cos;
    using boost::multiprecision::sin;
    vector<pair<BigFloat, BigFloat>> pts(n);
    const BigFloat pi = acos(BigFloat(-1));
    for (int i = 0; i < n; ++i) {
        BigFloat ang = BigFloat(2) * pi * i / n;
        pts[i] = {cos(ang), sin(ang)};
    }

    unordered_map<PointKey, int, PointKeyHash> counts;
    counts.reserve(static_cast<size_t>(choose4(n) + 1));
    counts.max_load_factor(0.7f);

    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            for (int c = b + 1; c < n; ++c) {
                for (int d = c + 1; d < n; ++d) {
                    const auto& A = pts[a];
                    const auto& C = pts[c];
                    const auto& B = pts[b];
                    const auto& D = pts[d];

                    BigFloat x1 = A.first, y1 = A.second;
                    BigFloat x2 = C.first, y2 = C.second;
                    BigFloat x3 = B.first, y3 = B.second;
                    BigFloat x4 = D.first, y4 = D.second;

                    BigFloat den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
                    BigFloat det1 = x1 * y2 - y1 * x2;
                    BigFloat det2 = x3 * y4 - y3 * x4;
                    BigFloat px = (det1 * (x3 - x4) - (x1 - x2) * det2) / den;
                    BigFloat py = (det1 * (y3 - y4) - (y1 - y2) * det2) / den;

                    // Quantize high-precision coordinates for stable hashing.
                    PointKey key{
                        round_to_int64(px),
                        round_to_int64(py)
                    };
                    ++counts[key];
                }
            }
        }
    }

    vector<int64_t> Nm(n / 2 + 2, 0);
    for (const auto& kv : counts) {
        int q = kv.second;
        int64_t D = 1 + 8LL * q;
        int64_t s = static_cast<int64_t>(llroundl(sqrtl(static_cast<long double>(D))));
        while (s * s < D) ++s;
        while (s * s > D) --s;
        if (s * s != D) {
            cerr << "Validation failure: non-triangular count for n=" << n << '\n';
            continue;
        }
        int m = static_cast<int>((1 + s) / 2);
        if (m >= static_cast<int>(Nm.size())) Nm.resize(m + 1, 0);
        Nm[m] += 1;
    }
    return Nm;
}

int64_t compute_T(int n, const Precomp& p) {
    vector<int64_t> Nm = intersection_multiplicities(n);
    int64_t total = 0;
    for (int m = 2; m < static_cast<int>(Nm.size()); ++m) {
        if (Nm[m] == 0) continue;
        int64_t ways = count_at_least_two(n, m, p);
        total = (total + mod_mul(Nm[m] % MOD, ways)) % MOD;
    }
    return total;
}

}  // namespace

int main() {
    const int N_MIN = 3;
    const int N_MAX = 60;
    Precomp pre = precompute(N_MAX);

    vector<int64_t> Tvals(N_MAX + 1, 0);
    int thread_count = static_cast<int>(thread::hardware_concurrency());
    if (thread_count <= 0) thread_count = 1;
    vector<int64_t> partial(thread_count, 0);

    atomic<int> next_n(N_MIN);
    vector<thread> workers;
    workers.reserve(thread_count);
    for (int t = 0; t < thread_count; ++t) {
        workers.emplace_back([&, t]() {
            for (;;) {
                int n = next_n.fetch_add(1);
                if (n > N_MAX) break;
                int64_t Tn = compute_T(n, pre);
                Tvals[n] = Tn;
                partial[t] = (partial[t] + Tn) % MOD;
            }
        });
    }
    for (auto& th : workers) th.join();

    if (Tvals[5] != 20 || Tvals[8] != 14640) {
        cerr << "Validation failure: T(5) or T(8) mismatch\n";
        return 1;
    }

    int64_t answer = 0;
    for (int t = 0; t < thread_count; ++t) {
        answer = (answer + partial[t]) % MOD;
    }
    cout << answer << '\n';
    return 0;
}
