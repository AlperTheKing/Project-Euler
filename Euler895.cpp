#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

static constexpr long long MOD = 989898989LL;

static long long add_mod(long long a, long long b) {
    long long v = a + b;
    if (v >= MOD) v -= MOD;
    if (v < 0) v += MOD;
    return v;
}

static long long mul_mod(long long a, long long b) {
    return static_cast<long long>((__int128)a * b % MOD);
}

static long long pair_sum_leq(int M1, int M2, int S) {
    if (S < 2) return 0;
    if (M1 > M2) swap(M1, M2);
    int max_sum = M1 + M2;
    if (S >= max_sum) return 1LL * M1 * M2;
    if (S <= M1 + 1) {
        long long s = S - 1;
        return s * (s + 1) / 2;
    }
    if (S <= M2 + 1) {
        long long first = 1LL * M1 * (M1 + 1) / 2;
        long long extra = 1LL * (S - (M1 + 1)) * M1;
        return first + extra;
    }
    // Descending tail.
    int k = max_sum - S;
    long long total = 1LL * M1 * M2;
    return total - 1LL * k * (k + 1) / 2;
}

static long long count_two_pos_one_neg(int Mp1, int Mp2, int Mn, int B) {
    if (Mp1 <= 0 || Mp2 <= 0 || Mn <= 0) return 0;
    int low = B + 1;
    int high = B + Mn;
    long long count = pair_sum_leq(Mp1, Mp2, high) - pair_sum_leq(Mp1, Mp2, low - 1);
    return max(0LL, count);
}

struct Precomp {
    int maxL;
    // D2[d][carry_idx][u_idx]
    vector<array<array<long long, 3>, 3>> d2;
    // D3slice[len][c0_idx][carry_idx][u_idx], u in [-6,6].
    vector<array<array<array<long long, 13>, 5>, 2>> d3slice;
    // F[len][s1_idx][c2_idx][u2_idx][u_final_idx]
    // s1_idx: 0 for -1, 1 for +1.
    vector<array<array<array<array<long long, 3>, 3>, 2>, 2>> f;
};

static Precomp build_precomp(int m) {
    Precomp pc;
    pc.maxL = max(0, m - 1);
    int maxL = pc.maxL;

    // Precompute D2.
    pc.d2.assign(maxL + 1, {});
    auto &d2 = pc.d2;
    // Indices: carry -1,0,1 -> 0,1,2. U -2,0,2 -> 0,1,2.
    d2[0][1][1] = 1; // carry=0, U=0.
    const int carry_vals[3] = {-1, 0, 1};
    const int u_vals[3] = {-2, 0, 2};
    const int s_vals[3] = {-2, 0, 2};
    const int s_mult[3] = {1, 2, 1};
    for (int d = 1; d <= maxL; ++d) {
        auto &cur = d2[d];
        auto &prev = d2[d - 1];
        for (int ci = 0; ci < 3; ++ci) {
            for (int ui = 0; ui < 3; ++ui) {
                long long cnt = prev[ci][ui];
                if (!cnt) continue;
                int c = carry_vals[ci];
                int u = u_vals[ui];
                for (int si = 0; si < 3; ++si) {
                    int S = s_vals[si];
                    int total = c + S;
                    if (total & 1) continue;
                    int c2 = total / 2;
                    int u2 = u + S;
                    int c2i = c2 + 1;
                    int u2i = (u2 + 2) / 2; // maps -2,0,2 -> 0,1,2.
                    cur[c2i][u2i] = add_mod(cur[c2i][u2i], mul_mod(cnt, s_mult[si]));
                }
            }
        }
    }

    // Precompute D3 slices.
    pc.d3slice.assign(maxL, {}); // lengths 0..maxL-1
    const int s3_vals[4] = {-3, -1, 1, 3};
    const int s3_mult[4] = {1, 3, 3, 1};
    for (int c0_idx = 0; c0_idx < 2; ++c0_idx) {
        int c0 = (c0_idx == 0) ? -1 : 1;
        int offset = 0;
        vector<array<long long, 5>> dist(1); // u range size 1 at l=0
        for (auto &row : dist) row.fill(0);
        dist[0][c0 + 2] = 1; // u=0 at offset 0.

        for (int len = 0; len < maxL; ++len) {
            // Store slice for this len.
            auto &slice = pc.d3slice[len][c0_idx];
            for (int ci = 0; ci < 5; ++ci) {
                for (int ui = 0; ui < 13; ++ui) {
                    int u = ui - 6;
                    if (abs(u) <= 3 * len) {
                        int idx = u + offset;
                        slice[ci][ui] = dist[idx][ci];
                    } else {
                        slice[ci][ui] = 0;
                    }
                }
            }

            if (len + 1 >= maxL) break;
            int new_offset = 3 * (len + 1);
            int new_size = 2 * new_offset + 1;
            vector<array<long long, 5>> next(static_cast<size_t>(new_size));
            for (auto &row : next) row.fill(0);

            int cur_size = 2 * offset + 1;
            for (int idx = 0; idx < cur_size; ++idx) {
                int u = idx - offset;
                for (int ci = 0; ci < 5; ++ci) {
                    long long cnt = dist[idx][ci];
                    if (!cnt) continue;
                    int c = ci - 2;
                    for (int si = 0; si < 4; ++si) {
                        int S = s3_vals[si];
                        int total = c + S;
                        if (total & 1) continue;
                        int c2 = total / 2;
                        int u2 = u + S;
                        int c2i = c2 + 2;
                        int idx2 = u2 + new_offset;
                        next[idx2][c2i] = add_mod(next[idx2][c2i], mul_mod(cnt, s3_mult[si]));
                    }
                }
            }
            dist.swap(next);
            offset = new_offset;
        }
    }

    // Precompute F for l >= 1 and S1 in {-1, +1}.
    pc.f.assign(maxL + 1, {});
    const int u2_vals[3] = {-2, 0, 2};
    const int u_final_vals[3] = {-1, 0, 1};
    for (int l = 1; l <= maxL; ++l) {
        for (int s1_idx = 0; s1_idx < 2; ++s1_idx) {
            int S1 = (s1_idx == 0) ? -1 : 1;
            auto &dst = pc.f[l][s1_idx];
            for (int c2_idx = 0; c2_idx < 2; ++c2_idx) {
                int c2 = (c2_idx == 0) ? -1 : 1;
                const auto &slice = pc.d3slice[l - 1][c2_idx];
                for (int u2_idx = 0; u2_idx < 3; ++u2_idx) {
                    int u2 = u2_vals[u2_idx];
                    for (int uf_idx = 0; uf_idx < 3; ++uf_idx) {
                        dst[c2_idx][u2_idx][uf_idx] = 0;
                    }
                    for (int c3_idx = 0; c3_idx < 5; ++c3_idx) {
                        int c3 = c3_idx - 2;
                        int total = c3 + S1;
                        if (total & 1) continue;
                        int carry_final = total / 2;
                        for (int u3_idx = 0; u3_idx < 13; ++u3_idx) {
                            long long cnt3 = slice[c3_idx][u3_idx];
                            if (!cnt3) continue;
                            int u3 = u3_idx - 6;
                            int U_pre = u2 + u3 + S1;
                            if (U_pre < -1 || U_pre > 1) continue;
                            if (carry_final != U_pre) continue;
                            int uf_idx = U_pre + 1; // -1->0,0->1,1->2
                            dst[c2_idx][u2_idx][uf_idx] = add_mod(dst[c2_idx][u2_idx][uf_idx], cnt3);
                        }
                    }
                }
            }
        }
    }

    return pc;
}

static long long compute_G(int m) {
    if (m <= 1) return 0;
    Precomp pc = build_precomp(m);
    int maxL = pc.maxL;

    long long total = 0;

    // Special case: all top lengths are 0.
    {
        int M = m;
        long long base = count_two_pos_one_neg(M, M, M, 0);
        long long add = mul_mod(6 % MOD, base % MOD); // 6 sign patterns with two positives / two negatives.
        total = add_mod(total, add);
    }

    const int u_final_vals[3] = {-1, 0, 1};
    const int carry_vals[3] = {-1, 0, 1};
    const int u2_vals[3] = {-2, 0, 2};

    auto worker = [&](int L_start, int L_end) -> long long {
        long long partial = 0;
        for (int L = L_start; L < L_end; ++L) {
            int M_long = m - L;
            for (int l = 0; l <= L; ++l) {
                int M_short = m - l;
                int d = L - l;

                long long top_neg[3] = {0, 0, 0};
                long long top_pos[3] = {0, 0, 0};

                if (l == 0) {
                    // Two-active only (long stacks). Compute for S1 = -2, 0, +2.
                    long long top_m2[3] = {0, 0, 0};
                    long long top_0[3] = {0, 0, 0};
                    long long top_p2[3] = {0, 0, 0};
                    const auto &d2 = pc.d2[L - 1];
                    for (int ci = 0; ci < 3; ++ci) {
                        for (int ui = 0; ui < 3; ++ui) {
                            long long cnt = d2[ci][ui];
                            if (!cnt) continue;
                            int c = carry_vals[ci];
                            int u = u2_vals[ui];
                            // S1 = -2
                            {
                                int total = c - 2;
                                if ((total & 1) == 0) {
                                    int carry_f = total / 2;
                                    int U_f = u - 2;
                                    if (U_f >= -1 && U_f <= 1 && carry_f == U_f) {
                                        top_m2[U_f + 1] = add_mod(top_m2[U_f + 1], cnt);
                                    }
                                }
                            }
                            // S1 = 0
                            {
                                int total = c;
                                if ((total & 1) == 0) {
                                    int carry_f = total / 2;
                                    int U_f = u;
                                    if (U_f >= -1 && U_f <= 1 && carry_f == U_f) {
                                        top_0[U_f + 1] = add_mod(top_0[U_f + 1], cnt);
                                    }
                                }
                            }
                            // S1 = +2
                            {
                                int total = c + 2;
                                if ((total & 1) == 0) {
                                    int carry_f = total / 2;
                                    int U_f = u + 2;
                                    if (U_f >= -1 && U_f <= 1 && carry_f == U_f) {
                                        top_p2[U_f + 1] = add_mod(top_p2[U_f + 1], cnt);
                                    }
                                }
                            }
                        }
                    }

                    for (int uf_idx = 0; uf_idx < 3; ++uf_idx) {
                        int U_final = u_final_vals[uf_idx];
                        int B = -U_final;
                        long long bottomA = count_two_pos_one_neg(M_long, M_long, M_short, B);
                        long long bottomB = count_two_pos_one_neg(M_long, M_short, M_long, B);
                        long long bottomA_neg = count_two_pos_one_neg(M_long, M_long, M_short, -B);
                        long long bottomB_neg = count_two_pos_one_neg(M_long, M_short, M_long, -B);

                        long long contrib_pos = add_mod(mul_mod(top_m2[uf_idx], bottomA % MOD),
                                                       mul_mod(top_0[uf_idx], (2LL * bottomB) % MOD));
                        long long contrib_neg = add_mod(mul_mod(top_p2[uf_idx], bottomA_neg % MOD),
                                                       mul_mod(top_0[uf_idx], (2LL * bottomB_neg) % MOD));
                        long long contrib = add_mod(contrib_pos, contrib_neg);
                        // 3 choices for which stack is short.
                        partial = add_mod(partial, mul_mod(3, contrib));
                    }
                    continue;
                }

                // l >= 1: use precomputed F for S1 = -1 and +1.
                const auto &d2 = pc.d2[d];
                const auto &f_neg = pc.f[l][0];
                const auto &f_pos = pc.f[l][1];
                for (int c2_idx = 0; c2_idx < 3; ++c2_idx) {
                    int c2 = carry_vals[c2_idx];
                    if (c2 == 0) continue; // no 3-active sequences from even carry.
                    int c2s_idx = (c2 == -1) ? 0 : 1;
                    for (int u2_idx = 0; u2_idx < 3; ++u2_idx) {
                        long long cnt2 = d2[c2_idx][u2_idx];
                        if (!cnt2) continue;
                        for (int uf_idx = 0; uf_idx < 3; ++uf_idx) {
                            long long add_neg = f_neg[c2s_idx][u2_idx][uf_idx];
                            if (add_neg) {
                                top_neg[uf_idx] = add_mod(top_neg[uf_idx], mul_mod(cnt2, add_neg));
                            }
                            long long add_pos = f_pos[c2s_idx][u2_idx][uf_idx];
                            if (add_pos) {
                                top_pos[uf_idx] = add_mod(top_pos[uf_idx], mul_mod(cnt2, add_pos));
                            }
                        }
                    }
                }

                for (int uf_idx = 0; uf_idx < 3; ++uf_idx) {
                    int U_final = u_final_vals[uf_idx];
                    int B = -U_final;
                    long long bottomA = count_two_pos_one_neg(M_long, M_long, M_short, B);
                    long long bottomB = count_two_pos_one_neg(M_long, M_short, M_long, B);
                    long long bottomA_neg = count_two_pos_one_neg(M_long, M_long, M_short, -B);
                    long long bottomB_neg = count_two_pos_one_neg(M_long, M_short, M_long, -B);

                    long long contrib_pos = add_mod(mul_mod(top_neg[uf_idx], (bottomA + 2LL * bottomB) % MOD), 0);
                    long long contrib_neg = add_mod(mul_mod(top_pos[uf_idx], (bottomA_neg + 2LL * bottomB_neg) % MOD), 0);
                    long long contrib = add_mod(contrib_pos, contrib_neg);
                    partial = add_mod(partial, mul_mod(3, contrib));
                }
            }
        }
        return partial;
    };

    int threads = thread::hardware_concurrency();
    if (threads <= 0) threads = 4;
    threads = min(threads, maxL);
    vector<long long> partials(threads, 0);
    vector<thread> pool;
    int chunk = (maxL + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
        int L_start = 1 + t * chunk;
        int L_end = min(maxL + 1, L_start + chunk);
        if (L_start >= L_end) continue;
        pool.emplace_back([&, t, L_start, L_end]() {
            partials[t] = worker(L_start, L_end);
        });
    }
    for (auto &th : pool) th.join();
    for (long long v : partials) {
        total = add_mod(total, v);
    }

    return total;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Validation checkpoints.
    struct Check {
        int m;
        long long expected;
    } checks[] = {
        {2, 6},
        {5, 348},
        {20, 125825982708LL}
    };
    bool ok = true;
    for (const auto &check : checks) {
        long long got = compute_G(check.m);
        long long exp = check.expected % MOD;
        if (got != exp) {
            cerr << "Validation failed for G(" << check.m << "): got " << got
                 << ", expected " << exp << ".\n";
            ok = false;
        }
    }
    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    } else {
        return 1;
    }

    const int target = 9898;
    long long answer = compute_G(target);
    cout << answer % MOD << "\n";
    return 0;
}
