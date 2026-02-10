#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

// Project Euler 572: Idempotent Matrices
//
// An integer 3x3 matrix M is idempotent iff it is a projection of the free
// abelian group Z^3, hence its rank over Q is 0,1,2,3. Rank 0 and 3 give
// the zero and identity matrices.
//
// Rank 1 idempotents have the form M = u v^T with u,v in Z^3 and v·u = 1.
// Rank 2 idempotents have the form M = I - u v^T with the same constraint.
// The pair (u,v) and (-u,-v) give the same matrix, so we count only vectors u
// with a canonical sign (first nonzero component positive) and count matching v.
//
// For rank 1, the bound |M_ij| <= n is equivalent to |u_i v_j| <= n for all i,j.
// For rank 2, bounds are |u_i v_j| <= n for i != j and |1 - u_i v_i| <= n.

using i64 = std::int64_t;
using u64 = std::uint64_t;

static inline int iabs(int x) { return x < 0 ? -x : x; }

static inline int gcd3(int a, int b, int c) {
    a = iabs(a);
    b = iabs(b);
    c = iabs(c);
    return std::gcd(a, std::gcd(b, c));
}

static inline bool canonical_u(int a, int b, int c) {
    if (a != 0) return a > 0;
    if (b != 0) return b > 0;
    return c > 0;  // (0,0,0) excluded by caller
}

static inline i64 floor_div(i64 a, i64 b) {
    // floor(a/b) for b != 0
    i64 q = a / b;
    i64 r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) --q;
    return q;
}

static inline i64 ceil_div(i64 a, i64 b) {
    // ceil(a/b) for b != 0
    return -floor_div(-a, b);
}

static inline u64 count_solutions_symmetric(int a, int b, int c, int V) {
    // Count integer (x,y,z) with -V<=x,y,z<=V and a*x + b*y + c*z = 1.
    if (V < 0) return 0;
    const int len = 2 * V + 1;

    const bool za = (a == 0), zb = (b == 0), zc = (c == 0);
    const int zcnt = int(za) + int(zb) + int(zc);

    if (zcnt == 2) {
        // One variable only.
        const int coef = !za ? a : (!zb ? b : c);
        if (coef == 0) return 0;
        if (1 % coef != 0) return 0;
        const int x = 1 / coef;
        if (iabs(x) > V) return 0;
        return static_cast<u64>(len) * static_cast<u64>(len);
    }
    if (zcnt == 1) {
        // Two variables.
        if (c == 0) {
            // a*x + b*y = 1, z free.
            u64 cnt = 0;
            for (int x = -V; x <= V; ++x) {
                const int rem = 1 - a * x;
                if (rem % b != 0) continue;
                const int y = rem / b;
                if (iabs(y) <= V) cnt += static_cast<u64>(len);
            }
            return cnt;
        }
        if (b == 0) {
            u64 cnt = 0;
            for (int x = -V; x <= V; ++x) {
                const int rem = 1 - a * x;
                if (rem % c != 0) continue;
                const int z = rem / c;
                if (iabs(z) <= V) cnt += static_cast<u64>(len);
            }
            return cnt;
        }
        // a==0
        u64 cnt = 0;
        for (int y = -V; y <= V; ++y) {
            const int rem = 1 - b * y;
            if (rem % c != 0) continue;
            const int z = rem / c;
            if (iabs(z) <= V) cnt += static_cast<u64>(len);
        }
        return cnt;
    }

    // Three-variable case. Pick pivot to minimize divisor checks: prefer +/-1, else largest |coef|.
    int coef[3] = {a, b, c};
    int pivot = 0;
    for (int k = 0; k < 3; ++k) {
        if (iabs(coef[k]) == 1) {
            pivot = k;
            break;
        }
        if (iabs(coef[k]) > iabs(coef[pivot])) pivot = k;
    }
    const int i = (pivot + 1) % 3;
    const int j = (pivot + 2) % 3;

    u64 cnt = 0;
    for (int xi = -V; xi <= V; ++xi) {
        for (int xj = -V; xj <= V; ++xj) {
            const int rem = 1 - coef[i] * xi - coef[j] * xj;
            const int den = coef[pivot];
            if (rem % den != 0) continue;
            const int xp = rem / den;
            if (iabs(xp) <= V) ++cnt;
        }
    }
    return cnt;
}

struct Range {
    int lo = 0;
    int hi = -1;
    int len() const { return hi >= lo ? (hi - lo + 1) : 0; }
};

static inline u64 count_solutions_box(int a, int b, int c, Range r0, Range r1, Range r2) {
    // Count (x,y,z) in r0 x r1 x r2 with a*x + b*y + c*z = 1.
    if (r0.len() == 0 || r1.len() == 0 || r2.len() == 0) return 0;

    int coef[3] = {a, b, c};
    Range rr[3] = {r0, r1, r2};

    // Handle degenerate cases with zero coefficients.
    const bool za = (a == 0), zb = (b == 0), zc = (c == 0);
    const int zcnt = int(za) + int(zb) + int(zc);
    if (zcnt == 3) return 0;

    if (zcnt == 2) {
        // One variable only.
        int k = !za ? 0 : (!zb ? 1 : 2);
        const int den = coef[k];
        if (1 % den != 0) return 0;
        const int x = 1 / den;
        return (x < rr[k].lo || x > rr[k].hi) ? 0ULL
                                              : static_cast<u64>(rr[(k + 1) % 3].len()) *
                                                    static_cast<u64>(rr[(k + 2) % 3].len());
    }

    if (zcnt == 1) {
        // Two variables + one free.
        int free_k = (za ? 0 : (zb ? 1 : 2));
        int i = (free_k + 1) % 3;
        int j = (free_k + 2) % 3;
        // Solve coef[i]*xi + coef[j]*xj = 1, x_free arbitrary.
        u64 cnt2 = 0;
        for (int xi = rr[i].lo; xi <= rr[i].hi; ++xi) {
            const int rem = 1 - coef[i] * xi;
            const int den = coef[j];
            if (rem % den != 0) continue;
            const int xj = rem / den;
            if (xj >= rr[j].lo && xj <= rr[j].hi) ++cnt2;
        }
        return cnt2 * static_cast<u64>(rr[free_k].len());
    }

    // Choose pivot with nonzero coefficient to minimize iterations on the other two ranges.
    int best_pivot = -1;
    i64 best_cost = (i64)1 << 62;
    for (int k = 0; k < 3; ++k) {
        if (coef[k] == 0) continue;
        int i = (k + 1) % 3;
        int j = (k + 2) % 3;
        i64 cost = static_cast<i64>(rr[i].len()) * static_cast<i64>(rr[j].len());
        // Prefer +/-1 pivots when tied.
        if (cost < best_cost || (cost == best_cost && iabs(coef[k]) == 1)) {
            best_cost = cost;
            best_pivot = k;
        }
    }
    const int k = best_pivot;
    const int i = (k + 1) % 3;
    const int j = (k + 2) % 3;

    u64 cnt = 0;
    for (int xi = rr[i].lo; xi <= rr[i].hi; ++xi) {
        for (int xj = rr[j].lo; xj <= rr[j].hi; ++xj) {
            const int rem = 1 - coef[i] * xi - coef[j] * xj;
            const int den = coef[k];
            if (rem % den != 0) continue;
            const int xk = rem / den;
            if (xk >= rr[k].lo && xk <= rr[k].hi) ++cnt;
        }
    }
    return cnt;
}

static inline Range intersect(Range a, Range b) {
    Range r;
    r.lo = std::max(a.lo, b.lo);
    r.hi = std::min(a.hi, b.hi);
    return r;
}

static inline u64 compute_C(int n) {
    // base: zero matrix always; identity only if n>=1
    const u64 base = 1 + (n >= 1 ? 1 : 0);

    const int threads = std::max(1u, std::min(8u, std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 4u));
    std::vector<std::thread> pool;
    std::vector<u64> part_r1(threads, 0), part_r2(threads, 0);

    const int a_min = -n;
    const int a_max = n;
    const int total_a = a_max - a_min + 1;
    const int block = (total_a + threads - 1) / threads;

    for (int t = 0; t < threads; ++t) {
        const int alo = a_min + t * block;
        const int ahi = std::min(a_max, alo + block - 1);
        pool.emplace_back([&, t, alo, ahi] {
            u64 r1 = 0, r2 = 0;
            for (int a = alo; a <= ahi; ++a) {
                for (int b = -n; b <= n; ++b) {
                    for (int c = -n; c <= n; ++c) {
                        if (a == 0 && b == 0 && c == 0) continue;
                        if (!canonical_u(a, b, c)) continue;
                        if (gcd3(a, b, c) != 1) continue;

                        const int Umax = std::max({iabs(a), iabs(b), iabs(c)});
                        const int V = n / Umax;
                        if (V > 0) {
                            r1 += count_solutions_symmetric(a, b, c, V);
                        }

                        // Rank 2: bounds for v0,v1,v2.
                        Range rv[3] = {{-n - 1, n + 1}, {-n - 1, n + 1}, {-n - 1, n + 1}};
                        const int u[3] = {a, b, c};
                        for (int j = 0; j < 3; ++j) {
                            int max_other = 0;
                            for (int i = 0; i < 3; ++i) {
                                if (i == j) continue;
                                max_other = std::max(max_other, iabs(u[i]));
                            }
                            if (max_other > 0) {
                                const int off = n / max_other;
                                rv[j] = intersect(rv[j], Range{-off, off});
                            }
                            if (u[j] != 0) {
                                // Diagonal constraint: 1-n <= u[j]*v[j] <= 1+n (note sign flip if u[j] < 0).
                                i64 lo, hi;
                                if (u[j] > 0) {
                                    lo = ceil_div(1LL - n, u[j]);
                                    hi = floor_div(1LL + n, u[j]);
                                } else {
                                    lo = ceil_div(1LL + n, u[j]);
                                    hi = floor_div(1LL - n, u[j]);
                                }
                                rv[j] = intersect(rv[j], Range{static_cast<int>(lo), static_cast<int>(hi)});
                            }
                        }
                        r2 += count_solutions_box(a, b, c, rv[0], rv[1], rv[2]);
                    }
                }
            }
            part_r1[t] = r1;
            part_r2[t] = r2;
        });
    }
    for (auto& th : pool) th.join();

    const u64 rank1 = std::accumulate(part_r1.begin(), part_r1.end(), 0ULL);
    const u64 rank2 = std::accumulate(part_r2.begin(), part_r2.end(), 0ULL);
    return base + rank1 + rank2;
}

int main() {
    // Validation from the statement.
    if (compute_C(1) != 164ULL) {
        std::cerr << "Validation failed: C(1)\n";
        return 1;
    }
    if (compute_C(2) != 848ULL) {
        std::cerr << "Validation failed: C(2)\n";
        return 1;
    }

    std::cout << compute_C(200) << "\n";
    return 0;
}
