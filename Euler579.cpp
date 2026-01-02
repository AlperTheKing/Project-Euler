#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;

namespace {

struct ValList {
    vector<int> v;
    vector<int> v2;
};

struct ValueSets {
    ValList even;
    ValList odd;
    ValList even_nonneg;
    ValList odd_nonneg;
};

struct Pattern {
    const ValList* a_vals;
    const ValList* b_vals;
    const ValList* c_vals;
    const ValList* d_vals;
    int g;
    int nmax;
};

struct Task {
    int pattern_idx;
    int a_start;
    int a_end;
};

inline int64_t iabs(int64_t x) {
    return x < 0 ? -x : x;
}

inline int gcd3(int64_t a, int64_t b, int64_t c) {
    return std::gcd((int)iabs(a), std::gcd((int)iabs(b), (int)iabs(c)));
}

ValueSets build_values(int M) {
    ValueSets out;
    out.even.v.reserve(M + 1);
    out.odd.v.reserve(M + 1);
    out.even_nonneg.v.reserve(M + 1);
    out.odd_nonneg.v.reserve(M + 1);

    for (int k = 0; k <= M; ++k) {
        if ((k & 1) == 0) {
            if (k == 0) {
                out.even.v.push_back(0);
                out.even_nonneg.v.push_back(0);
            } else {
                out.even.v.push_back(k);
                out.even.v.push_back(-k);
                out.even_nonneg.v.push_back(k);
            }
        } else {
            out.odd.v.push_back(k);
            out.odd.v.push_back(-k);
            out.odd_nonneg.v.push_back(k);
        }
    }

    auto fill_sq = [](ValList& list) {
        list.v2.resize(list.v.size());
        for (size_t i = 0; i < list.v.size(); ++i) {
            int x = list.v[i];
            list.v2[i] = x * x;
        }
    };

    fill_sq(out.even);
    fill_sq(out.odd);
    fill_sq(out.even_nonneg);
    fill_sq(out.odd_nonneg);
    return out;
}

unsigned __int128 compute_sum_frames(int n, unsigned int threads) {
    int m1 = (int)sqrt((long double)n);
    int m2 = (int)sqrt((long double)(2LL * n));
    int m3 = (int)sqrt((long double)(4LL * n));

    ValueSets vals1 = build_values(m1);
    ValueSets vals2 = build_values(m2);
    ValueSets vals3 = build_values(m3);

    vector<Pattern> patterns;
    patterns.reserve(15);

    for (int mask = 1; mask < 16; ++mask) {
        int odd_count = __builtin_popcount((unsigned)mask);
        if (odd_count == 0) continue;

        // For primitive quaternions, the frame gcd is 1/2/4 based on odd-component count.
        int g = 1;
        if (odd_count == 2) g = 2;
        else if (odd_count == 4) g = 4;

        int nmax = n * g;

        ValueSets* vs = nullptr;
        if (g == 1) vs = &vals1;
        else if (g == 2) vs = &vals2;
        else vs = &vals3;

        auto pick_full = [&](int bit) -> const ValList* {
            return (bit ? &vs->odd : &vs->even);
        };
        auto pick_nonneg = [&](int bit) -> const ValList* {
            return (bit ? &vs->odd_nonneg : &vs->even_nonneg);
        };

        Pattern pat;
        pat.a_vals = pick_nonneg(mask & 1);
        pat.b_vals = pick_full(mask & 2);
        pat.c_vals = pick_full(mask & 4);
        pat.d_vals = pick_full(mask & 8);
        pat.g = g;
        pat.nmax = nmax;
        patterns.push_back(pat);
    }

    vector<Task> tasks;
    const int chunk = 4;
    for (int idx = 0; idx < (int)patterns.size(); ++idx) {
        int a_size = (int)patterns[idx].a_vals->v.size();
        for (int i = 0; i < a_size; i += chunk) {
            Task t{idx, i, min(i + chunk, a_size)};
            tasks.push_back(t);
        }
    }

    if (threads == 0) threads = 1;
    vector<unsigned __int128> sums(threads, 0);
    atomic<size_t> task_index{0};

    auto worker = [&](unsigned int tid) {
        unsigned __int128 local_sum = 0;
        while (true) {
            size_t idx = task_index.fetch_add(1, memory_order_relaxed);
            if (idx >= tasks.size()) break;

            const Task& task = tasks[idx];
            const Pattern& pat = patterns[task.pattern_idx];
            const ValList& a_list = *pat.a_vals;
            const ValList& b_list = *pat.b_vals;
            const ValList& c_list = *pat.c_vals;
            const ValList& d_list = *pat.d_vals;

            int g = pat.g;
            int nmax = pat.nmax;

            for (int ai = task.a_start; ai < task.a_end; ++ai) {
                int a = a_list.v[ai];
                int a2 = a_list.v2[ai];
                bool a_zero = (a == 0);

                for (size_t bi = 0; bi < b_list.v.size(); ++bi) {
                    int b = b_list.v[bi];
                    if (a_zero && b < 0) continue;
                    int b2 = b_list.v2[bi];
                    int s1 = a2 + b2;
                    if (s1 > nmax) break;

                    int gcd_ab = std::gcd(iabs(a), iabs(b));
                    bool ab_zero = a_zero && (b == 0);

                    for (size_t ci = 0; ci < c_list.v.size(); ++ci) {
                        int c = c_list.v[ci];
                        if (ab_zero && c < 0) continue;
                        int c2 = c_list.v2[ci];
                        int s2 = s1 + c2;
                        if (s2 > nmax) break;

                        int gcd_abc = std::gcd(gcd_ab, (int)iabs(c));
                        bool abc_zero = ab_zero && (c == 0);

                        for (size_t di = 0; di < d_list.v.size(); ++di) {
                            int d = d_list.v[di];
                            if (abc_zero && d <= 0) continue;
                            int d2 = d_list.v2[di];
                            int N = s2 + d2;
                            if (N > nmax) break;
                            if (N == 0) continue;

                            if (std::gcd(gcd_abc, (int)iabs(d)) != 1) continue;

                            int64_t u1 = (int64_t)(a2 + b2 - c2 - d2) / g;
                            int64_t u2 = (int64_t)(2LL * (b * c - a * d)) / g;
                            int64_t u3 = (int64_t)(2LL * (b * d + a * c)) / g;

                            int64_t v1 = (int64_t)(2LL * (b * c + a * d)) / g;
                            int64_t v2 = (int64_t)(a2 - b2 + c2 - d2) / g;
                            int64_t v3 = (int64_t)(2LL * (c * d - a * b)) / g;

                            int64_t w1 = (int64_t)(2LL * (b * d - a * c)) / g;
                            int64_t w2 = (int64_t)(2LL * (c * d + a * b)) / g;
                            int64_t w3 = (int64_t)(a2 - b2 - c2 + d2) / g;

                            // Axis span equals the sum of absolute components for that axis.
                            int wx = (int)(iabs(u1) + iabs(v1) + iabs(w1));
                            int wy = (int)(iabs(u2) + iabs(v2) + iabs(w2));
                            int wz = (int)(iabs(u3) + iabs(v3) + iabs(w3));

                            if (wx > n || wy > n || wz > n) continue;
                            int maxw = max(wx, max(wy, wz));
                            int smax = n / maxw;
                            if (smax <= 0) continue;

                            int64_t L0 = (int64_t)N / g;
                            int gu = gcd3(u1, u2, u3);
                            int gv = gcd3(v1, v2, v3);
                            int gw = gcd3(w1, w2, w3);
                            int64_t G0 = (int64_t)gu + gv + gw;

                            int64_t ax = (int64_t)n + 1 - wx;
                            int64_t ay = (int64_t)n + 1 - wy;
                            int64_t az = (int64_t)n + 1 - wz;
                            int64_t L = L0;
                            int64_t G = G0;

                            for (int s = 1; s <= smax; ++s) {
                                unsigned long long placements = (unsigned long long)ax * (unsigned long long)ay *
                                                               (unsigned long long)az;

                                unsigned long long L64 = (unsigned long long)L;
                                unsigned long long G64 = (unsigned long long)G;
                                // Lattice points in a lattice cube: L^3 + (L+1) * (gu+gv+gw) + 1.
                                unsigned long long points = L64 * L64 * L64 + (L64 + 1ULL) * G64 + 1ULL;

                                local_sum += (unsigned __int128)placements * (unsigned __int128)points;

                                ax -= wx;
                                ay -= wy;
                                az -= wz;
                                L += L0;
                                G += G0;
                            }
                        }
                    }
                }
            }
        }
        sums[tid] = local_sum;
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (unsigned int t = 0; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for (auto& th : pool) th.join();

    unsigned __int128 total = 0;
    for (unsigned int t = 0; t < threads; ++t) total += sums[t];
    return total;
}

unsigned __int128 compute_S_exact(int n, unsigned int threads) {
    unsigned __int128 sum_frames = compute_sum_frames(n, threads);
    if (sum_frames % 24 != 0) {
        cerr << "Internal error: sum not divisible by 24.\n";
        exit(1);
    }
    return sum_frames / 24;
}

}  // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    auto check = [&](int n, unsigned long long expected) {
        unsigned __int128 got = compute_S_exact(n, 1);
        unsigned long long got64 = (unsigned long long)got;
        if (got64 != expected) {
            cerr << "Validation failed: S(" << n << ") got " << got64
                 << " expected " << expected << "\n";
            exit(1);
        }
    };

    check(1, 8ULL);
    check(2, 91ULL);
    check(4, 1878ULL);
    check(5, 5832ULL);
    check(10, 387003ULL);
    check(50, 29948928129ULL);
    cerr << "Validation checkpoints passed.\n";

    const int n = 5000;
    unsigned int threads = thread::hardware_concurrency();
    if (threads == 0) threads = 4;

    unsigned __int128 total = compute_S_exact(n, threads);
    unsigned long long mod = (unsigned long long)(total % 1000000000ULL);
    cout << mod << "\n";
    return 0;
}
