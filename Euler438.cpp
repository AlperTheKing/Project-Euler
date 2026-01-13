#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using int64 = long long;
using i128 = __int128_t;
using boost::multiprecision::cpp_int;

struct Interval {
    long double l;
    long double u;
    bool open_l;
    bool open_u;
};

struct Bounds {
    int64 lo;
    int64 hi;
};

static int64 factorial(int n) {
    int64 f = 1;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

static int64 comb(int n, int k) {
    if (k < 0 || k > n) return 0;
    int64 res = 1;
    for (int i = 1; i <= k; ++i) {
        res = res * (n - k + i) / i;
    }
    return res;
}

static bool is_integer_ld(long double v) {
    long double r = std::round(v);
    return fabsl(v - r) < 1e-12L;
}

static Bounds bounds_for_k(int n, int k, const std::vector<Interval>& roots) {
    struct Factor {
        long double lo;
        long double hi;
        bool open_lo;
        bool open_hi;
    };

    std::vector<Factor> factors;
    factors.reserve(n);
    for (const auto& r : roots) {
        long double lo = static_cast<long double>(k) - r.u;
        long double hi = static_cast<long double>(k) - r.l;
        bool open_lo = r.open_u;
        bool open_hi = r.open_l;
        factors.push_back({lo, hi, open_lo, open_hi});
    }

    long double min_val = 0.0L;
    long double max_val = 0.0L;
    bool min_att = false;
    bool max_att = false;
    bool first = true;

    int total = 1 << n;
    for (int mask = 0; mask < total; ++mask) {
        long double prod = 1.0L;
        bool att = true;
        for (int i = 0; i < n; ++i) {
            const auto& f = factors[i];
            if (mask & (1 << i)) {
                prod *= f.hi;
                if (f.open_hi) att = false;
            } else {
                prod *= f.lo;
                if (f.open_lo) att = false;
            }
        }
        if (first) {
            min_val = max_val = prod;
            min_att = max_att = att;
            first = false;
        } else {
            if (prod < min_val - 1e-18L) {
                min_val = prod;
                min_att = att;
            } else if (fabsl(prod - min_val) < 1e-18L && att) {
                min_att = true;
            }
            if (prod > max_val + 1e-18L) {
                max_val = prod;
                max_att = att;
            } else if (fabsl(prod - max_val) < 1e-18L && att) {
                max_att = true;
            }
        }
    }

    if (!min_att && is_integer_ld(min_val)) min_val += 1.0L;
    if (!max_att && is_integer_ld(max_val)) max_val -= 1.0L;

    int64 lo = static_cast<int64>(std::ceil(min_val - 1e-12L));
    int64 hi = static_cast<int64>(std::floor(max_val + 1e-12L));
    return {lo, hi};
}

static std::vector<Bounds> compute_bounds(int n, const std::vector<Interval>& roots) {
    std::vector<Bounds> bounds;
    bounds.reserve(n + 1);
    for (int k = 1; k <= n + 1; ++k) bounds.push_back(bounds_for_k(n, k, roots));
    return bounds;
}

static i128 eval_poly_int(const std::vector<int64>& coeffs, int64 x) {
    i128 res = 0;
    for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i) {
        res = res * x + coeffs[i];
    }
    return res;
}

static std::vector<cpp_int> poly_trim_cpp(std::vector<cpp_int> p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
    if (p.empty()) p.push_back(0);
    return p;
}

static std::vector<cpp_int> poly_derivative_cpp(const std::vector<cpp_int>& p) {
    if (p.size() <= 1) return {0};
    std::vector<cpp_int> d(p.size() - 1, 0);
    for (size_t i = 1; i < p.size(); ++i) d[i - 1] = p[i] * static_cast<long long>(i);
    return poly_trim_cpp(d);
}

static std::vector<cpp_int> poly_pseudo_remainder_cpp(std::vector<cpp_int> a, std::vector<cpp_int> b) {
    a = poly_trim_cpp(std::move(a));
    b = poly_trim_cpp(std::move(b));
    if (b.size() <= 1) return {0};
    if (b.back() < 0) {
        for (auto& v : b) v = -v;
    }

    int deg_a = static_cast<int>(a.size()) - 1;
    int deg_b = static_cast<int>(b.size()) - 1;
    if (deg_a < deg_b) return a;

    cpp_int lc = b.back();
    while (deg_a >= deg_b) {
        cpp_int lc_a = a.back();
        int shift = deg_a - deg_b;
        for (auto& v : a) v *= lc;
        for (int j = 0; j <= deg_b; ++j) a[j + shift] -= lc_a * b[j];
        while (!a.empty() && a.back() == 0) a.pop_back();
        deg_a = static_cast<int>(a.size()) - 1;
    }
    if (a.empty()) a.push_back(0);
    return a;
}

static std::vector<std::vector<cpp_int>> build_sturm_cpp(const std::vector<int64>& coeffs) {
    std::vector<cpp_int> p(coeffs.size(), 0);
    for (size_t i = 0; i < coeffs.size(); ++i) p[i] = static_cast<long long>(coeffs[i]);
    auto p0 = poly_trim_cpp(std::move(p));
    auto p1 = poly_derivative_cpp(p0);
    std::vector<std::vector<cpp_int>> seq;
    seq.push_back(p0);
    seq.push_back(p1);
    while (seq.back().size() > 1) {
        auto rem = poly_pseudo_remainder_cpp(seq[seq.size() - 2], seq.back());
        for (auto& v : rem) v = -v;
        seq.push_back(poly_trim_cpp(std::move(rem)));
    }
    return seq;
}

static cpp_int eval_poly_cpp(const std::vector<cpp_int>& p, int64 x) {
    cpp_int res = 0;
    for (int i = static_cast<int>(p.size()) - 1; i >= 0; --i) res = res * x + p[i];
    return res;
}

static int sign_variations_cpp(const std::vector<std::vector<cpp_int>>& sturm, int64 x) {
    int prev = 0;
    int changes = 0;
    for (const auto& poly : sturm) {
        cpp_int val = eval_poly_cpp(poly, x);
        int s = (val > 0) - (val < 0);
        if (s == 0) continue;
        if (prev != 0 && s != prev) ++changes;
        prev = s;
    }
    return changes;
}

static int count_interval(const std::vector<std::vector<cpp_int>>& sturm,
                          const std::vector<int64>& coeffs,
                          int64 a, int64 b,
                          bool include_left, bool include_right) {
    int va = sign_variations_cpp(sturm, a);
    int vb = sign_variations_cpp(sturm, b);
    int cnt = va - vb; // roots in (a, b]
    if (include_left && eval_poly_int(coeffs, a) == 0) cnt += 1;
    if (!include_right && eval_poly_int(coeffs, b) == 0) cnt -= 1;
    return cnt;
}

static bool check_intervals(int n, const std::vector<int64>& coeffs) {
    auto sturm = build_sturm_cpp(coeffs);
    for (int i = 1; i <= n; ++i) {
        int cnt = count_interval(sturm, coeffs, i, i + 1, true, false);
        if (cnt != 1) return false;
    }
    return true;
}

static bool check_prefix_divisibility(const std::vector<int64>& y,
                                      int idx,
                                      const std::vector<int64>& factorials) {
    std::vector<int64> diffs(y.begin(), y.begin() + idx + 1);
    for (int k = 1; k <= idx; ++k) {
        for (int i = 0; i <= idx - k; ++i) diffs[i] = diffs[i + 1] - diffs[i];
        if (diffs[0] % factorials[k] != 0) return false;
    }
    return true;
}

static bool build_coeffs(const std::vector<int64>& y,
                         const std::vector<int64>& factorials,
                         std::vector<int64>& coeffs,
                         i128& sum_abs) {
    int n = static_cast<int>(y.size()) - 1;
    std::vector<int64> diffs = y;
    std::vector<int64> d(n + 1, 0);

    for (int k = 0; k <= n; ++k) {
        if (k > 0) {
            for (int i = 0; i <= n - k; ++i) diffs[i] = diffs[i + 1] - diffs[i];
        }
        int64 f = factorials[k];
        if (diffs[0] % f != 0) return false;
        d[k] = diffs[0] / f;
    }

    std::vector<i128> poly(n + 1, 0);
    std::vector<i128> term(n + 1, 0);
    term[0] = 1;

    for (int k = 0; k <= n; ++k) {
        if (k > 0) {
            std::vector<i128> next(n + 1, 0);
            for (int deg = 0; deg <= k - 1; ++deg) {
                next[deg] -= static_cast<i128>(k) * term[deg];
                next[deg + 1] += term[deg];
            }
            term.swap(next);
        }
        if (d[k] != 0) {
            for (int deg = 0; deg <= k; ++deg) poly[deg] += static_cast<i128>(d[k]) * term[deg];
        }
    }

    if (poly[n] != 1) return false;

    coeffs.assign(n + 1, 0);
    sum_abs = 0;
    for (int i = 0; i <= n; ++i) coeffs[i] = static_cast<int64>(poly[i]);
    for (int i = 0; i < n; ++i) sum_abs += std::llabs(static_cast<long long>(coeffs[i]));
    return true;
}

struct Result {
    i128 sum = 0;
    int64 count = 0;
};

struct Solver {
    int n;
    std::vector<Bounds> bounds;
    std::vector<int64> coeffs;
    std::vector<int64> factorials;
    int64 fact_n;

    Result solve_range(int64 start, int64 end, std::atomic<int64>* progress) const {
        Result res;
        std::vector<int64> y(n + 1, 0);
        std::vector<int64> poly_coeffs;

        std::function<void(int, int64)> dfs = [&](int idx, int64 current_sum) {
            if (idx == n) {
                int64 val = fact_n - current_sum;
                if (val < bounds[n].lo || val > bounds[n].hi) return;
                for (int prev = 0; prev < idx; ++prev) {
                    int dist = idx - prev;
                    if ((val - y[prev]) % dist != 0) return;
                }
                y[idx] = val;
                if (!check_prefix_divisibility(y, idx, factorials)) return;
                i128 sum_abs = 0;
                if (!build_coeffs(y, factorials, poly_coeffs, sum_abs)) return;
                if (!check_intervals(n, poly_coeffs)) return;
                res.sum += sum_abs;
                res.count += 1;
                return;
            }

            int64 target = fact_n - current_sum;
            int64 min_rem = 0;
            int64 max_rem = 0;
            for (int k = idx; k <= n; ++k) {
                int64 c = coeffs[k];
                int64 lb = bounds[k].lo;
                int64 ub = bounds[k].hi;
                if (c >= 0) {
                    min_rem += c * lb;
                    max_rem += c * ub;
                } else {
                    min_rem += c * ub;
                    max_rem += c * lb;
                }
            }
            if (target < min_rem || target > max_rem) return;

            int64 lo = bounds[idx].lo;
            int64 hi = bounds[idx].hi;
            if (idx == 0) {
                lo = std::max(lo, start);
                hi = std::min(hi, end);
                if (lo > hi) return;
            }

            for (int64 val = lo; val <= hi; ++val) {
                if (idx == 0 && progress) {
                    progress->fetch_add(1, std::memory_order_relaxed);
                }
                bool ok = true;
                for (int prev = 0; prev < idx; ++prev) {
                    int dist = idx - prev;
                    if ((val - y[prev]) % dist != 0) { ok = false; break; }
                }
                if (!ok) continue;
                y[idx] = val;
                if (!check_prefix_divisibility(y, idx, factorials)) continue;
                dfs(idx + 1, current_sum + coeffs[idx] * val);
            }
        };

        dfs(0, 0);
        return res;
    }
};

static Solver make_solver(int n) {
    Solver solver;
    solver.n = n;
    solver.factorials.resize(n + 1, 1);
    for (int i = 1; i <= n; ++i) solver.factorials[i] = solver.factorials[i - 1] * i;
    solver.fact_n = solver.factorials[n];

    std::vector<Interval> roots;
    roots.reserve(n);
    for (int i = 1; i <= n; ++i) {
        roots.push_back({static_cast<long double>(i), static_cast<long double>(i + 1), false, true});
    }
    solver.bounds = compute_bounds(n, roots);
    solver.coeffs.resize(n + 1);
    for (int j = 0; j <= n; ++j) {
        int64 c = comb(n, j);
        if ((n - j) % 2 != 0) c = -c;
        solver.coeffs[j] = c;
    }
    return solver;
}

static std::string to_string_i128(i128 v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    if (neg) v = -v;
    std::string s;
    while (v > 0) {
        int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static std::string format_duration(double seconds) {
    if (seconds < 0) seconds = 0;
    int64 total = static_cast<int64>(seconds + 0.5);
    int64 h = total / 3600;
    int64 m = (total % 3600) / 60;
    int64 s = total % 60;
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << h << ":"
        << std::setw(2) << m << ":" << std::setw(2) << s;
    return oss.str();
}

static Result solve_parallel(const Solver& solver, int threads) {
    Result total;
    std::mutex mtx;

    int64 start = solver.bounds[0].lo;
    int64 end = solver.bounds[0].hi;
    int64 range = end - start + 1;
    if (range <= 0) return total;

    std::atomic<int64> progress{0};
    std::atomic<bool> done{false};
    auto start_time = std::chrono::steady_clock::now();
    std::thread progress_thread([&]() {
        size_t last_len = 0;
        while (!done.load(std::memory_order_relaxed)) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start_time).count();
            int64 done_count = progress.load(std::memory_order_relaxed);
            double pct = (range > 0) ? (100.0 * done_count / range) : 100.0;
            double rate = (elapsed > 0.0) ? (done_count / elapsed) : 0.0;
            std::ostringstream oss;
            oss << "Progress: " << std::fixed << std::setprecision(2)
                << pct << "% (" << done_count << "/" << range << ")"
                << " rate " << std::setprecision(1) << rate << "/s";
            std::string line = oss.str();
            if (line.size() < last_len) line.append(last_len - line.size(), ' ');
            last_len = line.size();
            std::cerr << '\r' << line << std::flush;
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    });

    if (threads <= 0) threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 4;
    threads = std::min<int64>(threads, range);

    int64 chunk = (range + threads - 1) / threads;
    std::vector<std::thread> pool;

    for (int t = 0; t < threads; ++t) {
        int64 s = start + t * chunk;
        int64 e = std::min(end, s + chunk - 1);
        if (s > e) continue;
        pool.emplace_back([&solver, &mtx, &total, &progress, s, e]() {
            Result local = solver.solve_range(s, e, &progress);
            std::lock_guard<std::mutex> lock(mtx);
            total.sum += local.sum;
            total.count += local.count;
        });
    }

    for (auto& th : pool) th.join();
    done.store(true, std::memory_order_relaxed);
    if (progress_thread.joinable()) progress_thread.join();
    {
        std::ostringstream oss;
        oss << "Progress: 100.00% (" << range << "/" << range << ")"
            << " rate 0.0/s";
        std::string line = oss.str();
        std::cerr << '\r' << line << '\n';
    }
    return total;
}

static bool run_validation() {
    std::cout << "Running validation for n=4...\n";
    Solver solver = make_solver(4);
    Result res = solve_parallel(solver, 1);
    std::cout << "n=4 count = " << res.count << ", sum = " << to_string_i128(res.sum) << "\n";
    if (res.count != 12 || res.sum != 2087) {
        std::cout << "Validation failed.\n";
        return false;
    }
    std::cout << "Validation passed.\n";
    std::cout << "-----------------------------\n";
    return true;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (!run_validation()) return 1;

    int n = 7;
    Solver solver = make_solver(n);
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    Result res = solve_parallel(solver, threads);

    std::cout << "n=7 count = " << res.count << "\n";
    std::cout << "Answer = " << to_string_i128(res.sum) << "\n";
    return 0;
}
