#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

// The squares game reduces to placing intervals of length 1 or sqrt(2) without overlap.
struct Segment {
    long double l;
    long double r;
    int grundy;
};

struct MSpline {
    long double l;
    long double r;
    long double slope;
    long double c;
};

struct Event {
    double pos;
    double delta;
    bool operator<(const Event& other) const {
        return pos < other.pos;
    }
};

static constexpr long double kEps = 1e-12L;

static std::vector<long double> build_breakpoints(long double limit, long double sqrt2) {
    std::vector<long double> vals;
    int max_b = static_cast<int>(limit / sqrt2) + 2;
    vals.reserve(static_cast<size_t>(limit * limit / 2));
    for (int b = 0; b <= max_b; ++b) {
        long double base = sqrt2 * static_cast<long double>(b);
        if (base > limit + kEps) break;
        int a_max = static_cast<int>(std::floor(limit - base + kEps));
        for (int a = 0; a <= a_max; ++a) {
            vals.push_back(base + static_cast<long double>(a));
        }
    }
    std::sort(vals.begin(), vals.end());
    std::vector<long double> uniq;
    uniq.reserve(vals.size());
    for (long double v : vals) {
        if (uniq.empty() || v - uniq.back() > kEps) {
            uniq.push_back(v);
        }
    }
    if (uniq.empty() || limit - uniq.back() > kEps) {
        uniq.push_back(limit);
    }
    return uniq;
}

static int find_segment(const std::vector<Segment>& segs, long double x) {
    int lo = 0;
    int hi = static_cast<int>(segs.size());
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (segs[mid].r <= x + kEps) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo >= static_cast<int>(segs.size())) return static_cast<int>(segs.size()) - 1;
    return lo;
}

static void collect_xors(long double u, const std::vector<Segment>& segs,
                         std::vector<char>& reachable) {
    if (u < -kEps) return;
    if (u <= kEps) {
        reachable[0] = 1;
        return;
    }
    if (segs.empty()) return;
    int i = 0;
    int j = find_segment(segs, u);
    while (i < static_cast<int>(segs.size()) && j >= 0) {
        long double x_l = segs[i].l;
        long double x_r = std::min(segs[i].r, u);
        if (x_l >= u - kEps) break;
        long double y_l = segs[j].l;
        long double y_r = segs[j].r;
        long double left = std::max(x_l, u - y_r);
        long double right = std::min(x_r, u - y_l);
        if (right > left + kEps) {
            int val = segs[i].grundy ^ segs[j].grundy;
            if (val >= static_cast<int>(reachable.size())) {
                reachable.resize(static_cast<size_t>(val + 1), 0);
            }
            reachable[val] = 1;
        }
        if (x_r <= u - y_l + kEps) {
            ++i;
        } else {
            --j;
        }
    }
}

static std::vector<Segment> build_grundy_segments(long double limit, long double sqrt2) {
    std::vector<long double> points = build_breakpoints(limit, sqrt2);
    std::vector<Segment> segs;
    std::vector<char> reachable(128, 0);
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        long double l = points[i];
        long double r = points[i + 1];
        long double mid = (l + r) * 0.5L;
        std::fill(reachable.begin(), reachable.end(), 0);
        if (mid >= 1.0L - kEps) {
            collect_xors(mid - 1.0L, segs, reachable);
        }
        if (mid >= sqrt2 - kEps) {
            collect_xors(mid - sqrt2, segs, reachable);
        }
        int mex = 0;
        while (mex < static_cast<int>(reachable.size()) && reachable[mex]) ++mex;
        if (mex >= static_cast<int>(reachable.size())) {
            reachable.resize(static_cast<size_t>(mex + 1), 0);
        }
        if (!segs.empty() && std::abs(segs.back().r - l) <= kEps && segs.back().grundy == mex) {
            segs.back().r = r;
        } else {
            segs.push_back({l, r, mex});
        }
    }
    return segs;
}

static std::vector<MSpline> build_m_splines(const std::vector<Segment>& gsegs,
                                            long double limit) {
    // m(u) = measure of x where Grundy(x) == Grundy(u-x); build as a piecewise-linear spline.
    int max_g = 0;
    for (const auto& seg : gsegs) max_g = std::max(max_g, seg.grundy);
    std::vector<std::vector<Segment>> per_g(static_cast<size_t>(max_g + 1));
    for (const auto& seg : gsegs) {
        per_g[seg.grundy].push_back(seg);
    }

    std::vector<Event> events;
    for (const auto& lst : per_g) {
        size_t n = lst.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                long double weight = (i == j) ? 1.0L : 2.0L;
                long double A = lst[i].l;
                long double B = lst[i].r;
                long double C = lst[j].l;
                long double D = lst[j].r;
                long double p1 = A + C;
                long double p2 = A + D;
                long double p3 = B + C;
                long double p4 = B + D;
                if (p1 <= limit + kEps) events.push_back({static_cast<double>(p1), static_cast<double>(weight)});
                if (p2 <= limit + kEps) events.push_back({static_cast<double>(p2), static_cast<double>(-weight)});
                if (p3 <= limit + kEps) events.push_back({static_cast<double>(p3), static_cast<double>(-weight)});
                if (p4 <= limit + kEps) events.push_back({static_cast<double>(p4), static_cast<double>(weight)});
            }
        }
    }

    std::sort(events.begin(), events.end());
    std::vector<MSpline> splines;
    long double prev = 0.0L;
    long double slope = 0.0L;
    long double value = 0.0L;
    size_t idx = 0;
    while (idx < events.size()) {
        long double pos = static_cast<long double>(events[idx].pos);
        if (pos > limit + kEps) break;
        if (pos > prev + kEps) {
            long double c = value - slope * prev;
            splines.push_back({prev, pos, slope, c});
            value += slope * (pos - prev);
            prev = pos;
        }
        long double delta = 0.0L;
        while (idx < events.size() && std::abs(static_cast<long double>(events[idx].pos) - pos) <= kEps) {
            delta += static_cast<long double>(events[idx].delta);
            ++idx;
        }
        slope += delta;
    }
    if (prev < limit - kEps) {
        long double c = value - slope * prev;
        splines.push_back({prev, limit, slope, c});
    }
    if (splines.empty()) {
        splines.push_back({0.0L, limit, 0.0L, 0.0L});
    }
    return splines;
}

static int find_m_segment(const std::vector<MSpline>& splines, long double x) {
    int lo = 0;
    int hi = static_cast<int>(splines.size());
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (splines[mid].r <= x + kEps) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo >= static_cast<int>(splines.size())) return static_cast<int>(splines.size()) - 1;
    return lo;
}

static long double eval_m(const std::vector<MSpline>& splines, long double u) {
    if (u <= 0.0L) return 0.0L;
    int idx = find_m_segment(splines, u);
    long double val = splines[idx].slope * u + splines[idx].c;
    if (val < 0.0L && val > -1e-9L) val = 0.0L;
    return val;
}

static long double e_value(const std::vector<MSpline>& splines, long double sqrt2, long double L) {
    long double u1 = L - 1.0L;
    long double u2 = L - sqrt2;
    if (u1 <= 0.0L || u2 <= 0.0L) return 0.0L;
    long double m1 = eval_m(splines, u1);
    long double m2 = eval_m(splines, u2);
    return 0.5L * L * (m1 / u1 + m2 / u2);
}

static long double fprime(long double L, long double sqrt2,
                          long double s1, long double c1,
                          long double s2, long double c2) {
    long double t1 = L - 1.0L;
    long double t2 = L - sqrt2;
    return (s1 + s2) - c1 / (t1 * t1) - c2 * sqrt2 / (t2 * t2);
}

static long double find_root(long double a, long double b, long double sqrt2,
                             long double s1, long double c1,
                             long double s2, long double c2) {
    long double fa = fprime(a, sqrt2, s1, c1, s2, c2);
    long double fb = fprime(b, sqrt2, s1, c1, s2, c2);
    if (std::abs(fa) <= 1e-18L) return a;
    if (std::abs(fb) <= 1e-18L) return b;
    if (fa * fb > 0.0L) return std::numeric_limits<long double>::quiet_NaN();
    for (int it = 0; it < 80; ++it) {
        long double mid = (a + b) * 0.5L;
        long double fm = fprime(mid, sqrt2, s1, c1, s2, c2);
        if (fa * fm <= 0.0L) {
            b = mid;
            fb = fm;
        } else {
            a = mid;
            fa = fm;
        }
    }
    return (a + b) * 0.5L;
}

static long double max_e_in_range(const std::vector<MSpline>& splines,
                                  const std::vector<long double>& breakpoints,
                                  long double sqrt2,
                                  long double Lmin, long double Lmax,
                                  unsigned threads) {
    int start = static_cast<int>(std::lower_bound(breakpoints.begin(), breakpoints.end(), Lmin - kEps) - breakpoints.begin());
    int end = static_cast<int>(std::upper_bound(breakpoints.begin(), breakpoints.end(), Lmax + kEps) - breakpoints.begin());
    int last_interval = static_cast<int>(breakpoints.size()) - 1;
    if (start >= last_interval) return 0.0L;
    if (end > last_interval) end = last_interval;
    if (end <= start) return 0.0L;

    int total = end - start;
    if (threads == 0) threads = 1;
    if (threads > static_cast<unsigned>(total)) threads = static_cast<unsigned>(total);

    std::vector<long double> best_per_thread(threads, 0.0L);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        int chunk_start = start + (total * static_cast<int>(t)) / static_cast<int>(threads);
        int chunk_end = start + (total * static_cast<int>(t + 1)) / static_cast<int>(threads);
        pool.emplace_back([&, t, chunk_start, chunk_end]() {
            long double local_best = 0.0L;
            for (int i = chunk_start; i < chunk_end; ++i) {
                long double a = std::max(breakpoints[i], Lmin);
                long double b = std::min(breakpoints[i + 1], Lmax);
                if (b <= a + kEps) continue;

                long double mid = (a + b) * 0.5L;
                long double u1 = mid - 1.0L;
                long double u2 = mid - sqrt2;
                if (u1 <= 0.0L || u2 <= 0.0L) continue;
                int idx1 = find_m_segment(splines, u1);
                int idx2 = find_m_segment(splines, u2);
                long double s1 = splines[idx1].slope;
                long double c1 = splines[idx1].c;
                long double s2 = splines[idx2].slope;
                long double c2 = splines[idx2].c;

                long double ea = e_value(splines, sqrt2, a);
                long double eb = e_value(splines, sqrt2, b);
                local_best = std::max(local_best, std::max(ea, eb));

                long double points[5] = {a, a + (b - a) * 0.25L, mid, a + (b - a) * 0.75L, b};
                long double fp[5];
                for (int p = 0; p < 5; ++p) {
                    fp[p] = fprime(points[p], sqrt2, s1, c1, s2, c2);
                }
                for (int p = 0; p < 4; ++p) {
                    if (fp[p] == 0.0L) {
                        long double ev = e_value(splines, sqrt2, points[p]);
                        local_best = std::max(local_best, ev);
                        continue;
                    }
                    if (fp[p] * fp[p + 1] <= 0.0L) {
                        long double root = find_root(points[p], points[p + 1], sqrt2, s1, c1, s2, c2);
                        if (root == root) {
                            long double ev = e_value(splines, sqrt2, root);
                            local_best = std::max(local_best, ev);
                        }
                    }
                }
            }
            best_per_thread[t] = local_best;
        });
    }
    for (auto& th : pool) th.join();
    long double best = 0.0L;
    for (long double v : best_per_thread) best = std::max(best, v);
    return best;
}

int main(int argc, char** argv) {
    long double Lmin = 200.0L;
    long double Lmax = 500.0L;
    if (argc > 2) {
        Lmin = std::strtold(argv[1], nullptr);
        Lmax = std::strtold(argv[2], nullptr);
    }
    long double sqrt2 = std::sqrt(2.0L);

    long double grundy_limit = Lmax;
    auto gsegs = build_grundy_segments(grundy_limit, sqrt2);

    long double u_limit = Lmax - 1.0L;
    auto splines = build_m_splines(gsegs, u_limit);

    long double global_min = std::min(Lmin, 2.0L);
    std::vector<long double> breakpoints;
    breakpoints.reserve(splines.size() * 4 + 8);
    breakpoints.push_back(global_min);
    breakpoints.push_back(Lmax);
    for (const auto& sp : splines) {
        long double a = sp.l + 1.0L;
        long double b = sp.r + 1.0L;
        long double c = sp.l + sqrt2;
        long double d = sp.r + sqrt2;
        if (a <= Lmax + kEps) breakpoints.push_back(a);
        if (b <= Lmax + kEps) breakpoints.push_back(b);
        if (c <= Lmax + kEps) breakpoints.push_back(c);
        if (d <= Lmax + kEps) breakpoints.push_back(d);
    }
    std::sort(breakpoints.begin(), breakpoints.end());
    std::vector<long double> uniq;
    uniq.reserve(breakpoints.size());
    for (long double v : breakpoints) {
        if (v < global_min - kEps || v > Lmax + kEps) continue;
        if (uniq.empty() || v - uniq.back() > kEps) uniq.push_back(v);
    }
    if (uniq.empty() || global_min - uniq.front() < -kEps) uniq.insert(uniq.begin(), global_min);
    if (uniq.empty() || Lmax - uniq.back() > kEps) uniq.push_back(Lmax);

    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    long double ans = max_e_in_range(splines, uniq, sqrt2, Lmin, Lmax, threads);

    // Validation checks.
    long double e2 = e_value(splines, sqrt2, 2.0L);
    if (std::abs(e2 - 2.0L) > 1e-7L) {
        std::cerr << "Validation failed: e(2) got " << std::fixed << std::setprecision(12) << e2
                  << " expected 2\n";
        return 1;
    }
    long double e4 = e_value(splines, sqrt2, 4.0L);
    if (std::abs(e4 - 1.11974851L) > 2e-6L) {
        std::cerr << "Validation failed: e(4) got " << std::fixed << std::setprecision(12) << e4
                  << " expected 1.11974851\n";
        return 1;
    }
    long double f2_10 = max_e_in_range(splines, uniq, sqrt2, 2.0L, 10.0L, threads);
    if (std::abs(f2_10 - 2.61969775L) > 2e-6L) {
        std::cerr << "Validation failed: f(2,10) got " << std::fixed << std::setprecision(12) << f2_10
                  << " expected 2.61969775\n";
        return 1;
    }
    long double f10_20 = max_e_in_range(splines, uniq, sqrt2, 10.0L, 20.0L, threads);
    if (std::abs(f10_20 - 5.99374121L) > 2e-6L) {
        std::cerr << "Validation failed: f(10,20) got " << std::fixed << std::setprecision(12) << f10_20
                  << " expected 5.99374121\n";
        return 1;
    }

    std::cout.setf(std::ios::fixed);
    std::cout.precision(8);
    std::cout << static_cast<double>(ans) << "\n";
    return 0;
}
