#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr long double kPerimeterEps = 1e-15L;

struct TriangleResult {
    long double perimeter = -1.0L;
    long double circumradius = -1.0L;
    i64 ax = 0;
    i64 ay = 0;
    i64 bx = 0;
    i64 by = 0;
    i64 cx = 0;
    i64 cy = 0;
};

struct Edge {
    i64 a = 0;
    i64 b = 0;
    u64 s2 = 0;
};

i128 floor_div(i128 a, i128 b) {
    i128 q = a / b;
    i128 r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) --q;
    return q;
}

i128 ceil_div(i128 a, i128 b) {
    i128 q = a / b;
    i128 r = a % b;
    if (r != 0 && ((r > 0) == (b > 0))) ++q;
    return q;
}

i128 round_div_nearest(i128 a, i128 b) {
    // Rounds a/b to nearest integer, ties away from zero.
    if (a >= 0) {
        return (a + b / 2) / b;
    }
    return -((-a + b / 2) / b);
}

u64 isqrt_u64(u64 n) {
    long double x = std::sqrt(static_cast<long double>(n));
    u64 r = static_cast<u64>(x);
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

u64 ceil_sqrt_u64(u64 n) {
    u64 r = isqrt_u64(n);
    if (r * r < n) ++r;
    return r;
}

u128 isqrt_u128(u128 n) {
    if (n == 0) return 0;
    long double approx = std::sqrt(static_cast<long double>(n));
    u128 x = static_cast<u128>(approx);
    if (x == 0) x = 1;
    while ((x + 1) * (x + 1) <= n) ++x;
    while (x * x > n) --x;
    return x;
}

i64 extended_gcd(i64 a, i64 b, i64 &x, i64 &y) {
    // Returns g=gcd(a,b), and x,y such that a*x + b*y = g.
    i64 x0 = 1, y0 = 0;
    i64 x1 = 0, y1 = 1;
    i64 aa = a, bb = b;
    while (bb != 0) {
        i64 q = aa / bb;
        i64 t = aa - q * bb;
        aa = bb;
        bb = t;

        i64 nx = x0 - q * x1;
        i64 ny = y0 - q * y1;
        x0 = x1;
        y0 = y1;
        x1 = nx;
        y1 = ny;
    }
    x = x0;
    y = y0;
    return aa;
}

long double perimeter_upper_bound(u64 s2) {
    long double L = std::sqrt(static_cast<long double>(s2));
    return 2.0L * L + 1.0L / L;
}

bool better_triangle(const TriangleResult &cand, const TriangleResult &best) {
    if (cand.perimeter > best.perimeter + kPerimeterEps) return true;
    if (std::fabsl(cand.perimeter - best.perimeter) <= kPerimeterEps &&
        cand.circumradius > best.circumradius + 1e-18L) {
        return true;
    }
    return false;
}

class Euler562Solver {
public:
    Euler562Solver(i64 r, int threads)
        : r_(r), r2_(static_cast<u64>(r) * static_cast<u64>(r)),
          threads_(std::max(1, threads)) {
        precompute_ymax();
    }

    TriangleResult solve_exact() {
        TriangleResult best;
        u64 d_bound = 128;
        const u64 max_d = 4ULL * r2_;

        while (true) {
            if (d_bound > max_d) d_bound = max_d;
            TriangleResult cur = search_with_deficiency_bound(d_bound, best);
            if (better_triangle(cur, best)) best = cur;

            if (best.perimeter < 0.0L) {
                if (d_bound == max_d) break;
                d_bound = std::min<u64>(max_d, d_bound * 2ULL);
                continue;
            }

            const long double p = best.perimeter;
            const long double disc = p * p - 8.0L;
            const long double l_bound = (p + std::sqrt(disc)) / 4.0L;
            const long double d_req_ld =
                4.0L * static_cast<long double>(r2_) - l_bound * l_bound;
            u64 d_req = 0;
            if (d_req_ld > 0.0L) {
                d_req = static_cast<u64>(std::ceill(d_req_ld - 1e-18L));
            }

            if (d_bound >= d_req || d_bound == max_d) break;
            d_bound = std::max<u64>(std::min<u64>(max_d, d_bound * 2ULL), d_req);
        }

        return best;
    }

private:
    i64 r_;
    u64 r2_;
    int threads_;
    std::vector<int> y_max_;

    void precompute_ymax() {
        y_max_.assign(static_cast<std::size_t>(2 * r_ + 1), 0);
        for (i64 x = -r_; x <= r_; ++x) {
            u64 rem = r2_ - static_cast<u64>(x * x);
            y_max_[static_cast<std::size_t>(x + r_)] = static_cast<int>(isqrt_u64(rem));
        }
    }

    TriangleResult search_with_deficiency_bound(u64 d_bound,
                                                const TriangleResult &initial_best) const {
        std::vector<Edge> edges;
        build_candidate_edges(d_bound, initial_best.perimeter, edges);
        if (edges.empty()) return initial_best;

        std::sort(edges.begin(), edges.end(), [](const Edge &lhs, const Edge &rhs) {
            if (lhs.s2 != rhs.s2) return lhs.s2 > rhs.s2;
            if (lhs.a != rhs.a) return lhs.a < rhs.a;
            return lhs.b < rhs.b;
        });

        TriangleResult best = initial_best;

        const std::size_t batch_size = 2048;
        std::size_t i = 0;
        while (i < edges.size()) {
            if (best.perimeter >= 0.0L &&
                perimeter_upper_bound(edges[i].s2) <= best.perimeter + kPerimeterEps) {
                break;
            }

            std::size_t j = std::min(edges.size(), i + batch_size);
            const int use_threads =
                (threads_ <= 1 || (j - i) < 64) ? 1 : std::min<int>(threads_, static_cast<int>(j - i));

            if (use_threads == 1) {
                for (std::size_t idx = i; idx < j; ++idx) {
                    const Edge &e = edges[idx];
                    if (best.perimeter >= 0.0L &&
                        perimeter_upper_bound(e.s2) <= best.perimeter + kPerimeterEps) {
                        continue;
                    }
                    TriangleResult cand = evaluate_edge(e, best.perimeter);
                    if (better_triangle(cand, best)) best = cand;
                }
            } else {
                std::vector<TriangleResult> locals(static_cast<std::size_t>(use_threads), best);
                std::vector<std::thread> pool;
                pool.reserve(static_cast<std::size_t>(use_threads));

                for (int t = 0; t < use_threads; ++t) {
                    pool.emplace_back([&, t]() {
                        TriangleResult local_best = best;
                        for (std::size_t idx = i + static_cast<std::size_t>(t); idx < j;
                             idx += static_cast<std::size_t>(use_threads)) {
                            const Edge &e = edges[idx];
                            if (local_best.perimeter >= 0.0L &&
                                perimeter_upper_bound(e.s2) <= local_best.perimeter + kPerimeterEps) {
                                continue;
                            }
                            TriangleResult cand = evaluate_edge(e, local_best.perimeter);
                            if (better_triangle(cand, local_best)) local_best = cand;
                        }
                        locals[static_cast<std::size_t>(t)] = local_best;
                    });
                }
                for (auto &th : pool) th.join();
                for (const auto &cand : locals) {
                    if (better_triangle(cand, best)) best = cand;
                }
            }

            i = j;
        }

        return best;
    }

    void build_candidate_edges(u64 d_bound, long double best_perimeter,
                               std::vector<Edge> &edges) const {
        edges.clear();
        const u64 max_s2 = 4ULL * r2_;
        const u64 min_s2 = (d_bound >= max_s2) ? 0ULL : (max_s2 - d_bound);

        // By symmetry we only need 1 <= b <= a.
        const u64 a_min_sq = (min_s2 + 1ULL) / 2ULL;
        i64 a_min = static_cast<i64>(ceil_sqrt_u64(a_min_sq));
        if (a_min < 1) a_min = 1;

        for (i64 a = a_min; a <= 2 * r_; ++a) {
            const u64 aa = static_cast<u64>(a) * static_cast<u64>(a);
            if (aa > max_s2) break;

            const u64 rem_hi = max_s2 - aa;
            i64 b_max = static_cast<i64>(isqrt_u64(rem_hi));
            if (b_max > a) b_max = a;
            if (b_max < 1) continue;

            i64 b_min = 1;
            if (aa < min_s2) {
                const u64 rem_lo = min_s2 - aa;
                b_min = static_cast<i64>(ceil_sqrt_u64(rem_lo));
            }
            if (b_min > b_max) continue;

            for (i64 b = b_min; b <= b_max; ++b) {
                if (std::gcd(a, b) != 1) continue;
                u64 s2 = aa + static_cast<u64>(b) * static_cast<u64>(b);
                if (s2 < min_s2 || s2 > max_s2) continue;
                if (best_perimeter >= 0.0L &&
                    perimeter_upper_bound(s2) <= best_perimeter + kPerimeterEps) {
                    continue;
                }
                edges.push_back({a, b, s2});
            }
        }
    }

    TriangleResult evaluate_edge(const Edge &e, long double threshold) const {
        TriangleResult best;
        best.perimeter = threshold;
        best.circumradius = -1.0L;

        const i64 a = e.a;
        const i64 b = e.b;
        const u64 s2_u64 = e.s2;
        const i128 s2 = static_cast<i128>(s2_u64);
        const long double L = std::sqrt(static_cast<long double>(s2_u64));

        if (perimeter_upper_bound(s2_u64) <= threshold + kPerimeterEps) return best;

        i64 x1 = 0, y1 = 0;
        i64 g = extended_gcd(-b, a, x1, y1);
        if (g == -1) {
            x1 = -x1;
            y1 = -y1;
        } else if (g != 1) {
            return best;
        }

        const i64 x_lo = std::max(-r_, -r_ - a);
        const i64 x_hi = std::min(r_, r_ - a);
        if (x_lo > x_hi) return best;

        for (i64 xA = x_lo; xA <= x_hi; ++xA) {
            const i64 xB = xA + a;
            const int y1_max = y_max_[static_cast<std::size_t>(xA + r_)];
            const int y2_max = y_max_[static_cast<std::size_t>(xB + r_)];
            const i64 y_lo = std::max<i64>(-y1_max, -static_cast<i64>(y2_max) - b);
            const i64 y_hi = std::min<i64>(y1_max, static_cast<i64>(y2_max) - b);
            if (y_lo > y_hi) continue;

            for (i64 yA = y_lo; yA <= y_hi; ++yA) {
                const i64 yB = yA + b;
                const i128 detA = static_cast<i128>(a) * yA - static_cast<i128>(b) * xA;

                for (int sgn : {-1, 1}) {
                    const i128 k = detA + static_cast<i128>(sgn);

                    // Raw point on line det(u, C)=k: C_raw = k * (x1, y1).
                    i128 x_raw = static_cast<i128>(x1) * k;
                    i128 y_raw = static_cast<i128>(y1) * k;

                    // Shift along u to nearest point to origin to keep numbers bounded.
                    i128 dot = x_raw * a + y_raw * b;
                    i128 shift = round_div_nearest(-dot, s2);
                    i128 x0 = x_raw + shift * a;
                    i128 y0 = y_raw + shift * b;

                    const i128 Acoef = s2;
                    const i128 Bcoef = 2 * (x0 * a + y0 * b);
                    const i128 Ccoef = x0 * x0 + y0 * y0 - static_cast<i128>(r2_);

                    i128 disc = Bcoef * Bcoef - 4 * Acoef * Ccoef;
                    if (disc < 0) continue;

                    u128 sqrt_disc = isqrt_u128(static_cast<u128>(disc));
                    i128 denom = 2 * Acoef;
                    i128 t_min = ceil_div(-Bcoef - static_cast<i128>(sqrt_disc), denom);
                    i128 t_max = floor_div(-Bcoef + static_cast<i128>(sqrt_disc), denom);
                    if (t_min > t_max) continue;

                    for (i128 t = t_min; t <= t_max; ++t) {
                        i128 xC128 = x0 + t * a;
                        i128 yC128 = y0 + t * b;
                        i128 normC2 = xC128 * xC128 + yC128 * yC128;
                        if (normC2 > static_cast<i128>(r2_)) continue;

                        i64 xC = static_cast<i64>(xC128);
                        i64 yC = static_cast<i64>(yC128);

                        if ((xC == xA && yC == yA) || (xC == xB && yC == yB)) continue;

                        i128 dx1 = xC128 - static_cast<i128>(xA);
                        i128 dy1 = yC128 - static_cast<i128>(yA);
                        i128 d2_128 = dx1 * dx1 + dy1 * dy1;

                        i128 dx2 = xC128 - static_cast<i128>(xB);
                        i128 dy2 = yC128 - static_cast<i128>(yB);
                        i128 d3_128 = dx2 * dx2 + dy2 * dy2;

                        if (d2_128 > s2 || d3_128 > s2) continue;  // u is longest side.

                        u64 d2 = static_cast<u64>(d2_128);
                        u64 d3 = static_cast<u64>(d3_128);

                        long double p =
                            L + std::sqrt(static_cast<long double>(d2)) + std::sqrt(static_cast<long double>(d3));
                        if (p <= best.perimeter + kPerimeterEps) continue;

                        long double R =
                            std::sqrt(static_cast<long double>(s2_u64) * static_cast<long double>(d2) *
                                      static_cast<long double>(d3)) /
                            2.0L;

                        TriangleResult cand;
                        cand.perimeter = p;
                        cand.circumradius = R;
                        cand.ax = xA;
                        cand.ay = yA;
                        cand.bx = xB;
                        cand.by = yB;
                        cand.cx = xC;
                        cand.cy = yC;

                        if (better_triangle(cand, best)) best = cand;
                    }
                }
            }
        }

        return best;
    }
};

long double compute_T(i64 r, int threads, TriangleResult *out_triangle = nullptr) {
    Euler562Solver solver(r, threads);
    TriangleResult best = solver.solve_exact();
    if (out_triangle) *out_triangle = best;
    return best.circumradius / static_cast<long double>(r);
}

bool run_validation(int threads) {
    struct Checkpoint {
        i64 r;
        long double expected;
        long double tol;
    };

    const std::vector<Checkpoint> checks = {
        {5, std::sqrt(19669.0L / 50.0L), 1e-12L},
        {10, 97.26729L, 1e-5L},
        {100, 9157.64707L, 1e-5L},
    };

    for (const auto &cp : checks) {
        TriangleResult tri;
        long double got = compute_T(cp.r, threads, &tri);
        long double err = std::fabsl(got - cp.expected);
        std::cout << "T(" << cp.r << ") = " << std::setprecision(12) << std::fixed << got << "\n";
        if (err > cp.tol) {
            std::cerr << "Validation failed at r=" << cp.r << ": expected ~" << cp.expected
                      << ", got " << got << "\n";
            std::cerr << "Triangle: A=(" << tri.ax << "," << tri.ay << "), B=(" << tri.bx << ","
                      << tri.by << "), C=(" << tri.cx << "," << tri.cy << ")\n";
            return false;
        }
    }

    std::cerr << "Validation checkpoints passed.\n";
    return true;
}

}  // namespace

int main(int argc, char **argv) {
    i64 target_r = 10'000'000;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--r" || arg == "-r") && i + 1 < argc) {
            target_r = std::stoll(argv[++i]);
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--no-validation") {
            validate = false;
        }
    }

    if (validate) {
        if (!run_validation(std::min(threads, 4))) return 1;
    }

    TriangleResult tri;
    long double T = compute_T(target_r, threads, &tri);
    unsigned long long rounded = static_cast<unsigned long long>(std::llround(T));

    std::cout << std::setprecision(12) << std::fixed;
    std::cout << "T(" << target_r << ") = " << T << "\n";
    std::cout << "Rounded answer = " << rounded << "\n";
    std::cout << rounded << "\n";

    return 0;
}
