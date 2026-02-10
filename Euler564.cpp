#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

namespace {

using i64 = std::int64_t;

static i64 round6(long double x) {
    return static_cast<i64>(std::llround(x * 1'000'000.0L));
}

struct Counts {
    // excess[r] = how many sides have length (r+1), for r>=1; length 1 is handled separately.
    std::vector<int> excess;
    int used_parts = 0;
    int max_len = 1;
};

struct EvalContext {
    int n = 0;
    int E = 0;
    long double log_total = 0;
    const std::vector<long double>* logfact = nullptr;
};

static void eval_partition(const Counts& cnt, const EvalContext& ctx, int c1, long double& expected_area) {
    const int n = ctx.n;
    const int max_len = cnt.max_len;
    const int perimeter = 2 * n - 3;

    // Compute log(weight) for this multiset: n! / (c1! * prod_r excess[r]!)
    long double logw = (*ctx.logfact)[n] - (*ctx.logfact)[c1];
    for (int r = 1; r <= ctx.E; ++r) {
        const int c = cnt.excess[r];
        if (c) logw -= (*ctx.logfact)[c];
    }
    const long double prob = expl(logw - ctx.log_total);

    const long double pi = acosl(-1.0L);

    // For a cyclic polygon of circumradius R and consecutive chord lengths l_i, define
    // alpha_i(R) = asin(l_i/(2R)) (half of the *minor* central angle).
    // If the circumcenter lies inside (or on) the polygon, all sides use minor arcs and:
    //   sum alpha_i(R) = pi.
    // If that equation has no solution, the maximal-area cyclic polygon has a unique longest
    // side that uses the major arc; then:
    //   alpha_max(R) = sum_{i!=max} alpha_i(R)  <=>  sum alpha_i(R) = 2*alpha_max(R).
    auto sums_at = [&](long double R) {
        const long double twoR = 2.0L * R;
        const long double fourR2 = 4.0L * R * R;
        long double sum_alpha = 0.0L;
        long double sum_dalpha = 0.0L;
        long double alpha_max = 0.0L;
        long double dalpha_max = 0.0L;

        // length 1 part (never the unique maximum when E>0, but included for completeness)
        if (c1) {
            const long double l = 1.0L;
            const long double s = sqrtl(fourR2 - l * l);
            const long double a = asinl(l / twoR);
            const long double da = -(l / (R * s));
            sum_alpha += static_cast<long double>(c1) * a;
            sum_dalpha += static_cast<long double>(c1) * da;
        }

        for (int r = 1; r <= ctx.E; ++r) {
            const int c = cnt.excess[r];
            if (!c) continue;
            const long double l = static_cast<long double>(r + 1);
            const long double s = sqrtl(fourR2 - l * l);
            const long double a = asinl(l / twoR);
            const long double da = -(l / (R * s));
            sum_alpha += static_cast<long double>(c) * a;
            sum_dalpha += static_cast<long double>(c) * da;
            if (r + 1 == max_len) {
                alpha_max = a;
                dalpha_max = da;
            }
        }

        return std::tuple<long double, long double, long double, long double>(sum_alpha, sum_dalpha, alpha_max, dalpha_max);
    };

    // Decide whether the "all minor arcs" solution exists by checking the minimal admissible radius.
    const long double low = 0.5L * static_cast<long double>(max_len);
    const auto [sum_alpha_low, sum_dalpha_low, alpha_max_low, dalpha_max_low] = sums_at(low);
    (void)sum_dalpha_low;
    (void)dalpha_max_low;

    const bool minor_exists = (sum_alpha_low >= pi - 1e-30L);
    const bool use_major = !minor_exists;
    if (use_major) {
        // With a major arc, the longest side must be unique; otherwise the major side cannot
        // have half-angle larger than the sum of the others.
        const int cmax = (max_len == 1) ? c1 : cnt.excess[max_len - 1];
        assert(cmax == 1);
        (void)alpha_max_low;  // silence unused in release builds
    }

    auto f_and_fp = [&](long double R) {
        const auto [sum_alpha, sum_dalpha, alpha_max, dalpha_max] = sums_at(R);
        if (!use_major) {
            return std::pair<long double, long double>(sum_alpha - pi, sum_dalpha);
        }
        return std::pair<long double, long double>(sum_alpha - 2.0L * alpha_max, sum_dalpha - 2.0L * dalpha_max);
    };

    long double R = low;
    if (!use_major) {
        // Root can sit exactly at low (e.g. when the longest side is a diameter).
        if (fabsl(sum_alpha_low - pi) > 1e-30L) {
            long double lo = low;
            long double hi = std::max(low * 2.0L, static_cast<long double>(perimeter) / (2.0L * pi));
            for (int it = 0; it < 200; ++it) {
                const auto [f, _fp] = f_and_fp(hi);
                if (f < 0) break;
                hi *= 2.0L;
            }

            R = 0.5L * (lo + hi);
            for (int it = 0; it < 40; ++it) {
                const auto [f, fp] = f_and_fp(R);
                if (fabsl(f) < 1e-24L) break;
                long double Rn = R - f / fp;
                if (!(Rn > lo && Rn < hi) || !std::isfinite(Rn)) Rn = 0.5L * (lo + hi);
                const auto [fn, _] = f_and_fp(Rn);
                if (fn > 0) {
                    // f(lo) > 0, f(hi) < 0
                    lo = Rn;
                } else {
                    hi = Rn;
                }
                R = Rn;
            }
        }
    } else {
        long double lo = low;
        long double hi = std::max(low * 2.0L, static_cast<long double>(perimeter) / (2.0L * pi));
        for (int it = 0; it < 200; ++it) {
            const auto [f, _fp] = f_and_fp(hi);
            if (f > 0) break;
            hi *= 2.0L;
        }
        R = 0.5L * (lo + hi);
        for (int it = 0; it < 60; ++it) {
            const auto [f, fp] = f_and_fp(R);
            if (fabsl(f) < 1e-24L) break;
            long double Rn = R - f / fp;
            if (!(Rn > lo && Rn < hi) || !std::isfinite(Rn)) Rn = 0.5L * (lo + hi);
            const auto [fn, _] = f_and_fp(Rn);
            if (fn < 0) {
                // f(lo) < 0, f(hi) > 0
                lo = Rn;
            } else {
                hi = Rn;
            }
            R = Rn;
        }
    }

    // Area of cyclic polygon: sum over sides of (l/4)*sqrt(4R^2 - l^2).
    const long double fourR2 = 4.0L * R * R;
    long double area_sum = 0.0L;
    long double contrib_max = 0.0L;
    if (c1) {
        const long double l = 1.0L;
        area_sum += static_cast<long double>(c1) * (l * 0.25L) * sqrtl(fourR2 - l * l);
    }
    for (int r = 1; r <= ctx.E; ++r) {
        const int c = cnt.excess[r];
        if (!c) continue;
        const long double l = static_cast<long double>(r + 1);
        const long double contrib = (l * 0.25L) * sqrtl(fourR2 - l * l);
        area_sum += static_cast<long double>(c) * contrib;
        if (r + 1 == max_len) contrib_max = contrib;
    }

    // If we used a major arc on the unique longest side, its contribution is negative, i.e.
    // area = sum_all - 2*contrib_max.
    const long double area = use_major ? fabsl(area_sum - 2.0L * contrib_max) : area_sum;
    expected_area += prob * area;
}

static long double E_of_n(int n, const std::vector<long double>& logfact) {
    const int E = n - 3;
    assert(E >= 0);

    // total compositions = C(E+n-1, n-1) = C(2n-4, n-1)
    const int m = 2 * n - 4;
    const long double log_total = logfact[m] - logfact[n - 1] - logfact[n - 3];

    long double expected_area = 0.0L;
    Counts cnt;
    cnt.excess.assign(E + 1, 0);
    cnt.used_parts = 0;
    cnt.max_len = 1;

    EvalContext ctx;
    ctx.n = n;
    ctx.E = E;
    ctx.log_total = log_total;
    ctx.logfact = &logfact;

    if (E == 0) {
        // All sides are 1: regular n-gon (here n=3) with side 1.
        eval_partition(cnt, ctx, n, expected_area);
        return expected_area;
    }

    auto rec = [&](auto&& self, int max_part, int remaining) -> void {
        if (remaining == 0) {
            eval_partition(cnt, ctx, n - cnt.used_parts, expected_area);
            return;
        }
        for (int part = std::min(max_part, remaining); part >= 1; --part) {
            cnt.excess[part] += 1;
            cnt.used_parts += 1;
            cnt.max_len = std::max(cnt.max_len, part + 1);
            self(self, part, remaining - part);
            cnt.excess[part] -= 1;
            cnt.used_parts -= 1;
            if (cnt.excess[part] == 0 && cnt.max_len == part + 1) {
                int ml = 1;
                for (int r = 1; r <= E; ++r) {
                    if (cnt.excess[r]) ml = std::max(ml, r + 1);
                }
                cnt.max_len = ml;
            }
        }
    };
    rec(rec, E, E);

    return expected_area;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Precompute log-factorials up to 100.
    std::vector<long double> logfact(101, 0.0L);
    for (int i = 1; i <= 100; ++i) logfact[i] = lgammal(static_cast<long double>(i) + 1.0L);

    long double S = 0.0L;
    for (int n = 3; n <= 50; ++n) {
        const long double En = E_of_n(n, logfact);
        S += En;

        // Statement validation points (rounded to 6 decimals).
        if (n == 3) assert(round6(En) == 433013);
        if (n == 4) assert(round6(En) == 1299038);
        if (n == 5) assert(round6(S) == 4604767);
        if (n == 10) assert(round6(S) == 66955511);
    }

    std::cout << std::fixed << std::setprecision(6) << static_cast<double>(S) << '\n';
    return 0;
}
