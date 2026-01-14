#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

static inline int depth_int(int m) {
    if (m <= 0) return -1;
    return 31 - __builtin_clz(static_cast<unsigned>(m));
}

struct BandRMQ {
    int d = 0;
    int start = 0;
    int end = 0;
    int len = 0;
    vector<long long> s1;
    vector<long long> s2;
    vector<int> lg;
    vector<vector<int>> st1;
    vector<vector<int>> st2;

    void build(int d_, int start_, int end_,
               const vector<long long>& A,
               const vector<long long>& B) {
        d = d_;
        start = start_;
        end = end_;
        len = end - start + 1;

        s1.resize(len);
        s2.resize(len);
        for (int i = 0; i < len; i++) {
            int m = start + i;
            s1[i] = A[m] - 1LL * m * (d + 1);
            s2[i] = B[m] - 1LL * m * d;
        }

        lg.assign(len + 1, 0);
        for (int i = 2; i <= len; i++) lg[i] = lg[i / 2] + 1;

        int K = lg[len] + 1;
        st1.assign(K, vector<int>(len));
        st2.assign(K, vector<int>(len));
        for (int i = 0; i < len; i++) {
            st1[0][i] = i;
            st2[0][i] = i;
        }
        for (int k = 1; k < K; k++) {
            int span = 1 << (k - 1);
            int limit = len - (1 << k) + 1;
            for (int i = 0; i < limit; i++) {
                int i1 = st1[k - 1][i];
                int i2 = st1[k - 1][i + span];
                st1[k][i] = (s1[i1] <= s1[i2]) ? i1 : i2;

                int j1 = st2[k - 1][i];
                int j2 = st2[k - 1][i + span];
                st2[k][i] = (s2[j1] <= s2[j2]) ? j1 : j2;
            }
        }
    }

    inline int rmqArgMin1Abs(int lAbs, int rAbs) const {
        int l = lAbs - start;
        int r = rAbs - start;
        int L = r - l + 1;
        int k = lg[L];
        int i1 = st1[k][l];
        int i2 = st1[k][r - (1 << k) + 1];
        int idx = (s1[i1] <= s1[i2]) ? i1 : i2;
        return start + idx;
    }

    inline int rmqArgMin2Abs(int lAbs, int rAbs) const {
        int l = lAbs - start;
        int r = rAbs - start;
        int L = r - l + 1;
        int k = lg[L];
        int i1 = st2[k][l];
        int i2 = st2[k][r - (1 << k) + 1];
        int idx = (s2[i1] <= s2[i2]) ? i1 : i2;
        return start + idx;
    }
};

static vector<long long> compute_exact_small(int n) {
    const long long INF = numeric_limits<long long>::max() / 4;
    vector<vector<long long>> dp(n + 2, vector<long long>(n + 2, 0));
    for (int len = 2; len <= n; ++len) {
        for (int l = 1; l + len - 1 <= n; ++l) {
            int r = l + len - 1;
            long long best = INF;
            for (int k = l; k <= r; ++k) {
                long long cost = static_cast<long long>(k) +
                                 max(dp[l][k - 1], dp[k + 1][r]);
                if (cost < best) best = cost;
            }
            dp[l][r] = best;
        }
    }
    vector<long long> out(n + 1, 0);
    for (int i = 1; i <= n; ++i) out[i] = dp[1][i];
    return out;
}

int main() {
    const int N = 200000;
    const long long INF = (1LL << 62);

    int maxD = depth_int(N);

    // ---- Compute A(m) ----
    vector<long long> A(N + 1, 0), B(N + 1, 0);

    vector<vector<long long>> sufMinVal(maxD + 1);
    vector<vector<int>> sufArgVal(maxD + 1);
    vector<char> built(maxD + 1, 0);
    built[0] = 1;

    auto buildSuf = [&](int d) {
        int lo = 1 << (d - 1);
        int hi = (1 << d) - 1;
        int len = hi - lo + 1;
        sufMinVal[d].assign(len, 0);
        sufArgVal[d].assign(len, 0);

        long long best = A[hi] - 1LL * hi * d;
        int arg = hi;
        sufMinVal[d][len - 1] = best;
        sufArgVal[d][len - 1] = arg;

        for (int i = len - 2; i >= 0; i--) {
            int q = lo + i;
            long long v = A[q] - 1LL * q * d;
            if (v < best) {
                best = v;
                arg = q;
            }
            sufMinVal[d][i] = best;
            sufArgVal[d][i] = arg;
        }
        built[d] = 1;
    };

    for (int m = 2; m <= N; m++) {
        int d = depth_int(m);
        if (d >= 1 && !built[d]) buildSuf(d);

        if (d == 0) {
            A[m] = 0;
            continue;
        }
        int base = 1 << (d - 1);
        int qLo = max(base, m - (1 << d));
        long long rightCost = 1LL * m * d + sufMinVal[d][qLo - base];

        long long leftCost = INF;
        int p = m - base + 1;
        if (p <= (1 << d)) {
            leftCost = 1LL * p + A[p - 1];
        }
        A[m] = min(leftCost, rightCost);
    }

    for (int d = 1; d <= maxD; d++) {
        if (!built[d]) buildSuf(d);
    }

    // ---- Compute B(m) ----
    auto B_cand_for_q = [&](int m, int q) -> long long {
        int d = depth_int(m);
        int p = m - q;
        int L = p - 1;
        long long best = 1LL * p * (d - 1) + B[q];
        int Dl = depth_int(L);
        if (Dl == d - 1) {
            best = max(best, 1LL * p + B[L]);
        } else if (Dl == d - 2) {
            best = max(best, 1LL * p + A[L]);
        }
        return best;
    };

    for (int m = 2; m <= N; m++) {
        int d = depth_int(m);
        if (d == 0) {
            B[m] = 0;
            continue;
        }
        int base = 1 << (d - 1);
        int p0 = m - base + 1;
        long long leftA = INF;
        if (p0 <= (1 << d)) {
            leftA = 1LL * p0 + A[p0 - 1];
        }

        int qLo = max(base, m - (1 << d));
        long long rightBestA = 1LL * m * d + sufMinVal[d][qLo - base];

        if (leftA == A[m]) {
            int L = p0 - 1;
            int R = m - p0;
            long long bestB = max(1LL * p0 + B[L], 1LL * p0 * (d - 1) + A[R]);

            if (rightBestA == A[m]) {
                int q = sufArgVal[d][qLo - base];
                long long cand = B_cand_for_q(m, q);
                bestB = min(bestB, cand);
            }
            B[m] = bestB;
        } else {
            int q = sufArgVal[d][qLo - base];
            B[m] = B_cand_for_q(m, q);
        }
    }

    // ---- Build RMQ per depth band ----
    vector<BandRMQ> bands(maxD + 1);
    vector<char> hasBand(maxD + 1, 0);
    for (int d = 1; d <= maxD; d++) {
        int start = 1 << d;
        int end = min((1 << (d + 1)) - 1, N);
        if (start <= end) {
            bands[d].build(d, start, end, A, B);
            hasBand[d] = 1;
        }
    }

    // ---- Compute C(n) and sum ----
    vector<long long> C(N + 1, 0);
    long long sum = 0;

    auto obj = [&](int n, int d, int m) -> long long {
        int k = n - m;
        long long L = 1LL * k + C[k - 1];
        long long g1 = 1LL * k * (d + 1) + A[m];
        long long g2 = 1LL * k * d + B[m];
        return max(L, max(g1, g2));
    };

    for (int n = 2; n <= N; n++) {
        long long best = INF;

        // d=0 band: m=1
        {
            int m = 1;
            int k = n - m;
            long long cand = 1LL * k + C[k - 1];
            best = min(best, cand);
        }

        int maxDn = depth_int(n - 1);
        for (int d = 1; d <= maxDn; d++) {
            int mLo = 1 << d;
            int mHi = min((1 << (d + 1)) - 1, n - 1);
            if (mLo > mHi) continue;

            auto rightCost = [&](int m) -> long long {
                int k = n - m;
                long long t1 = 1LL * k * d + A[m];
                long long t2 = 1LL * k * (d - 1) + B[m];
                return max(t1, t2);
            };
            auto P = [&](int m) -> bool {
                int k = n - m;
                return C[k - 1] <= rightCost(m);
            };
            auto h = [&](int m) -> long long {
                int k = n - m;
                return 1LL * k + A[m] - B[m];
            };

            if (!P(mHi)) {
                for (int m : {mHi, mHi - 1}) {
                    if (m < mLo || m > mHi) continue;
                    best = min(best, obj(n, d, m));
                }
                continue;
            }

            int mP;
            if (P(mLo)) {
                mP = mLo;
            } else {
                int lo = mLo, hi = mHi;
                while (lo < hi) {
                    int mid = (lo + hi) >> 1;
                    if (P(mid)) hi = mid; else lo = mid + 1;
                }
                mP = lo;
                for (int m : {mP - 1, mP}) {
                    if (m < mLo || m > mHi) continue;
                    best = min(best, obj(n, d, m));
                }
            }

            int subLo = mP;
            int subHi = mHi;

            for (int m : {subLo, subHi, subHi - 1, subHi - 2}) {
                if (m < subLo || m > subHi) continue;
                best = min(best, obj(n, d, m));
            }

            long long hLo = h(subLo);
            long long hHi = h(subHi);

            if (hasBand[d]) {
                const BandRMQ& band = bands[d];

                if (hLo >= 0 && hHi >= 0) {
                    int mMin = band.rmqArgMin1Abs(subLo, subHi);
                    for (int m = mMin - 2; m <= mMin + 2; m++) {
                        if (m < subLo || m > subHi) continue;
                        best = min(best, obj(n, d, m));
                    }
                } else if (hLo <= 0 && hHi <= 0) {
                    int mMin = band.rmqArgMin2Abs(subLo, subHi);
                    for (int m = mMin - 2; m <= mMin + 2; m++) {
                        if (m < subLo || m > subHi) continue;
                        best = min(best, obj(n, d, m));
                    }
                } else {
                    int lo = subLo, hi = subHi;
                    while (lo < hi) {
                        int mid = (lo + hi) >> 1;
                        if (h(mid) >= 0) hi = mid; else lo = mid + 1;
                    }
                    int mQ = lo;
                    for (int m = mQ - 4; m <= mQ + 4; m++) {
                        if (m < subLo || m > subHi) continue;
                        best = min(best, obj(n, d, m));
                    }

                    int mMin1 = band.rmqArgMin1Abs(subLo, subHi);
                    int mMin2 = band.rmqArgMin2Abs(subLo, subHi);
                    for (int mm : {mMin1, mMin2}) {
                        for (int m = mm - 2; m <= mm + 2; m++) {
                            if (m < subLo || m > subHi) continue;
                            best = min(best, obj(n, d, m));
                        }
                    }
                }
            }
        }

        C[n] = best;
        sum += best;
    }

    // ---- Validation checkpoints ----
    {
        auto exact = compute_exact_small(200);
        for (int i = 1; i <= 200; ++i) {
            assert(C[i] == exact[i]);
        }
    }
    assert(C[1] == 0);
    assert(C[2] == 1);
    assert(C[3] == 2);
    assert(C[8] == 12);
    assert(C[100] == 400);
    long long sum100 = 0;
    for (int i = 1; i <= 100; ++i) sum100 += C[i];
    assert(sum100 == 17575);
    assert(sum == 260511850222LL);

    cout << sum << "\n";
    return 0;
}
