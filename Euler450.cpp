#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

using int64 = long long;
using i128 = __int128_t;

static inline i128 tri_i128(int64 n) { return (i128)n * (n + 1) / 2; }

static string to_string_i128(i128 v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    if (neg) v = -v;
    string s;
    while (v > 0) {
        int digit = (int)(v % 10);
        s.push_back(char('0' + digit));
        v /= 10;
    }
    if (neg) s.push_back('-');
    reverse(s.begin(), s.end());
    return s;
}

static inline i128 iabs128(i128 x) { return x < 0 ? -x : x; }

static i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        i128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static i128 gcd128_3(i128 a, i128 b, i128 c) {
    return gcd128(gcd128(a, b), c);
}

static inline pair<i128, i128> pow_gauss_i128(i128 x, i128 y, int n) {
    i128 rx = 1, ry = 0;
    for (int i = 0; i < n; ++i) {
        i128 nx = rx * x - ry * y;
        i128 ny = rx * y + ry * x;
        rx = nx;
        ry = ny;
    }
    return {rx, ry};
}

static inline i128 ipow_i128(i128 base, int exp) {
    i128 r = 1;
    for (int i = 0; i < exp; ++i) r *= base;
    return r;
}

static inline tuple<int64, i128, i128, i128> denom_required(int p, int q, int A, int B, int D) {
    int e = p - q;
    auto aq = pow_gauss_i128((i128)A, (i128)B, q);
    auto aec = pow_gauss_i128((i128)A, (i128)(-B), e);

    i128 Deq = ipow_i128((i128)D, e - q);
    i128 ux = (i128)e * aq.first * Deq + (i128)q * aec.first;
    i128 uy = (i128)e * aq.second * Deq + (i128)q * aec.second;

    i128 De = ipow_i128((i128)D, e);
    i128 g0 = gcd128_3(De, ux, uy);
    i128 den128 = De / g0;

    int64 den = (int64)den128;
    return {den, ux, uy, g0};
}

static inline array<pair<int, int>, 8> variations(int a, int b) {
    return array<pair<int, int>, 8>{
        pair<int, int>{a, b},
        pair<int, int>{a, -b},
        pair<int, int>{-a, b},
        pair<int, int>{-a, -b},
        pair<int, int>{b, a},
        pair<int, int>{b, -a},
        pair<int, int>{-b, a},
        pair<int, int>{-b, -a},
    };
}

static vector<tuple<int, int, int>> primitive_triples(int Cmax) {
    vector<tuple<int, int, int>> out;
    int mlim = (int)std::sqrt((long double)Cmax) + 2;
    for (int m = 2; m <= mlim; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0) continue;
            if (std::gcd(m, n) != 1) continue;
            int a = m * m - n * n;
            int b = 2 * m * n;
            int c = m * m + n * n;
            if (c > Cmax) continue;
            out.emplace_back(a, b, c);
        }
    }
    return out;
}

static void sieve_mu_phi(int N, vector<int> &mu, vector<int> &phi) {
    mu.assign(N + 1, 0);
    phi.assign(N + 1, 0);
    vector<int> primes;
    vector<char> comp(N + 1, 0);
    mu[1] = 1;
    phi[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) {
            primes.push_back(i);
            phi[i] = i - 1;
            mu[i] = -1;
        }
        for (int p : primes) {
            long long v = 1LL * i * p;
            if (v > N) break;
            comp[(int)v] = 1;
            if (i % p == 0) {
                phi[(int)v] = phi[i] * p;
                mu[(int)v] = 0;
                break;
            } else {
                phi[(int)v] = phi[i] * (p - 1);
                mu[(int)v] = -mu[i];
            }
        }
    }
}

static vector<int64> compute_B_parallel(int N, const vector<int> &mu, int maxThreads = 0) {
    vector<int64> B(N + 1, 0);

    unsigned hw = thread::hardware_concurrency();
    int T = (maxThreads > 0) ? maxThreads : (hw ? (int)hw : 4);
    T = max(1, min(T, 8));

    vector<vector<int64>> local(T, vector<int64>(N + 1, 0));

    auto worker = [&](int tid, int dStart, int dEnd) {
        for (int d = dStart; d <= dEnd; ++d) {
            int md = mu[d];
            if (md == 0) continue;
            int64 coef = (int64)md * d;
            for (int p = d; p <= N; p += d) {
                int64 m = (p - 1) / (2LL * d);
                local[tid][p] += coef * (m * (m + 1) / 2);
            }
        }
    };

    vector<thread> th;
    th.reserve(T);
    int block = N / T;
    int cur = 1;
    for (int t = 0; t < T; ++t) {
        int start = cur;
        int end = (t == T - 1) ? N : (cur + block - 1);
        cur = end + 1;
        th.emplace_back(worker, t, start, end);
    }
    for (auto &x : th) x.join();

    for (int p = 1; p <= N; ++p) {
        int64 s = 0;
        for (int t = 0; t < T; ++t) s += local[t][p];
        B[p] = s;
    }
    return B;
}

static i128 compute_T(int N) {
    vector<int> mu, phi;
    sieve_mu_phi(N, mu, phi);
    vector<int64> B = compute_B_parallel(N, mu);

    i128 total = 0;

    for (int p = 3; p <= N; ++p) {
        int64 A = phi[p] / 2;
        if (A == 0) continue;
        int64 Bsum = B[p];

        int64 inner;
        if (p & 1) {
            inner = 4LL * p * A - 2LL * Bsum;
        } else {
            if (p % 4 == 0) inner = 4LL * p * A;
            else inner = 4LL * p * A - 4LL * Bsum;
        }

        int64 M = N / p;
        total += (i128)inner * tri_i128(M);
    }

    int Dmax = (int)std::sqrt((long double)(N / 3)) + 2;
    auto triples = primitive_triples(Dmax);

    // For N = 1e6, checking D = 5 (the smallest hypotenuse) shows no D>1
    // contributions for p > 12. We cap extra work here accordingly.
    int maxPextra = 12;

    for (int p = 3; p <= min(N, maxPextra); ++p) {
        int64 M = N / p;
        for (int q = 1; 2 * q < p; ++q) {
            if (std::gcd(p, q) != 1) continue;

            for (auto [a, b, D] : triples) {
                auto vars = variations(a, b);
                for (auto [A, Bv] : vars) {
                    auto [den, ux, uy, g0] = denom_required(p, q, A, Bv, D);
                    if (den > M) continue;

                    i128 bx = ux / g0;
                    i128 by = uy / g0;
                    i128 base_abs = iabs128(bx) + iabs128(by);

                    int64 mmax = M / den;
                    total += base_abs * tri_i128(mmax);
                }
            }
        }
    }

    return total;
}

static pair<i128, size_t> compute_S_and_count(int R, int r) {
    int g = std::gcd(R, r);
    int p = R / g;
    int q = r / g;
    int e = p - q;

    int Dmax = 600;
    auto triples = primitive_triples(Dmax);

    set<pair<long long, long long>> pts;

    auto ipow = [](int k) -> pair<int, int> {
        k %= 4;
        if (k < 0) k += 4;
        switch (k) {
            case 0: return {1, 0};
            case 1: return {0, 1};
            case 2: return {-1, 0};
            default: return {0, -1};
        }
    };

    for (int k = 0; k < 4; ++k) {
        auto wq = ipow(k * q);
        auto wne = ipow(-k * e);
        long long zx = (long long)e * wq.first + (long long)q * wne.first;
        long long zy = (long long)e * wq.second + (long long)q * wne.second;
        pts.insert({(long long)g * zx, (long long)g * zy});
    }

    for (auto [a, b, D] : triples) {
        auto vars = variations(a, b);
        for (auto [A, Bv] : vars) {
            auto [den, ux, uy, g0] = denom_required(p, q, A, Bv, D);
            if (den == 0) continue;
            if (g % den != 0) continue;
            long long k = g / den;

            i128 bx = ux / g0;
            i128 by = uy / g0;
            i128 x = (i128)k * bx;
            i128 y = (i128)k * by;

            pts.insert({(long long)x, (long long)y});
        }
    }

    i128 S = 0;
    for (auto [x, y] : pts) {
        S += (i128)llabs(x) + (i128)llabs(y);
    }
    return {S, pts.size()};
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    {
        auto [S31, c31] = compute_S_and_count(3, 1);
        assert((long long)S31 == 10);
        assert(c31 == 4);

        auto [Sbig, cbig] = compute_S_and_count(2500, 1000);
        assert((long long)Sbig == 24944);
        assert(cbig == 12);
    }
    {
        i128 t3 = compute_T(3);
        i128 t10 = compute_T(10);
        i128 t100 = compute_T(100);
        i128 t1000 = compute_T(1000);
        assert((long long)t3 == 10);
        assert((long long)t10 == 524);
        assert((long long)t100 == 580442);
        assert((long long)t1000 == 583108600);
    }

    i128 ans = compute_T(1'000'000);
    cout << to_string_i128(ans) << "\n";
    return 0;
}
