#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace std;

namespace {

const long double kPi = acosl(-1.0L);

long double clamp01(long double v) {
    if (v < 0.0L && v > -1e-15L) return 0.0L;
    if (v > 1.0L && v < 1.0L + 1e-15L) return 1.0L;
    return v;
}

long double H(int a, int b, long double x) {
    long double aa = static_cast<long double>(a);
    long double bb = static_cast<long double>(b);
    long double denom = 2.0L * (aa + bb);
    long double val = 0.5L - (bb * cosl(aa * kPi * x) + aa * cosl(bb * kPi * x)) / denom;
    return clamp01(val);
}

int derivative_sign(int a, int b, long double x) {
    int diff = abs(b - a);
    long double s = sinl((a + b) * kPi * x / 2.0L) * cosl(diff * kPi * x / 2.0L);
    if (s > 0.0L) return 1;
    if (s < 0.0L) return -1;
    return 0;
}

vector<long double> turning_values(int a, int b) {
    // Keep only points where the derivative flips sign (true extrema).
    int sum = a + b;
    int diff = abs(b - a);
    vector<long double> pts;
    pts.reserve(sum / 2 + diff / 2 + 2);

    pts.push_back(0.0L);
    pts.push_back(1.0L);

    for (int k = 1; k < sum / 2; ++k) {
        pts.push_back(2.0L * k / sum);
    }
    if (diff > 0) {
        for (int k = 0; k < diff / 2; ++k) {
            pts.push_back((2.0L * k + 1.0L) / diff);
        }
    }

    sort(pts.begin(), pts.end());
    vector<long double> uniq;
    uniq.reserve(pts.size());
    for (long double x : pts) {
        if (uniq.empty() || fabsl(x - uniq.back()) > 1e-18L) {
            uniq.push_back(x);
        }
    }

    vector<int> signs;
    signs.reserve(uniq.size());
    for (size_t i = 0; i + 1 < uniq.size(); ++i) {
        long double mid = (uniq[i] + uniq[i + 1]) / 2.0L;
        int s = derivative_sign(a, b, mid);
        if (s == 0) {
            mid = (uniq[i] * 0.75L + uniq[i + 1] * 0.25L);
            s = derivative_sign(a, b, mid);
        }
        signs.push_back(s);
    }

    vector<long double> turning;
    turning.reserve(uniq.size());
    turning.push_back(uniq.front());
    for (size_t i = 1; i + 1 < uniq.size(); ++i) {
        if (signs[i - 1] != signs[i]) {
            turning.push_back(uniq[i]);
        }
    }
    turning.push_back(uniq.back());

    vector<long double> values;
    values.reserve(turning.size());
    for (long double x : turning) {
        values.push_back(H(a, b, x));
    }
    return values;
}

long double path_variation(const vector<long double>& za, const vector<long double>& zb) {
    // Traverse the unique path by hopping between shared endpoints of monotone segments.
    int i = 0;
    int j = 0;
    long double zcur = 0.0L;
    int dir = 1;
    long double total = 0.0L;
    const long double eps = 1e-15L;
    long long max_steps = static_cast<long long>((za.size() + zb.size()) * (za.size() + zb.size()) + 10);
    long long steps = 0;

    while (steps < max_steps) {
        if (i == static_cast<int>(za.size()) - 1 && j == static_cast<int>(zb.size()) - 1) {
            return total;
        }
        if (i < 0 || j < 0 || i >= static_cast<int>(za.size()) - 1 || j >= static_cast<int>(zb.size()) - 1) {
            return numeric_limits<long double>::quiet_NaN();
        }

        long double fa0 = za[i];
        long double fa1 = za[i + 1];
        long double ga0 = zb[j];
        long double ga1 = zb[j + 1];
        long double min_f = min(fa0, fa1);
        long double max_f = max(fa0, fa1);
        long double min_g = min(ga0, ga1);
        long double max_g = max(ga0, ga1);

        long double next_z;
        bool hit_f = false;
        bool hit_g = false;

        if (dir > 0) {
            next_z = min(max_f, max_g);
            hit_f = fabsl(max_f - next_z) <= eps;
            hit_g = fabsl(max_g - next_z) <= eps;
            total += fabsl(next_z - zcur);
            zcur = next_z;
            if (hit_f) i += (fa1 >= fa0) ? 1 : -1;
            if (hit_g) j += (ga1 >= ga0) ? 1 : -1;
        } else {
            next_z = max(min_f, min_g);
            hit_f = fabsl(min_f - next_z) <= eps;
            hit_g = fabsl(min_g - next_z) <= eps;
            total += fabsl(next_z - zcur);
            zcur = next_z;
            if (hit_f) i += (fa1 <= fa0) ? 1 : -1;
            if (hit_g) j += (ga1 <= ga0) ? 1 : -1;
        }
        dir = -dir;
        ++steps;
    }

    return numeric_limits<long double>::quiet_NaN();
}

long double compute_F(int a, int b, int c, int d) {
    vector<long double> za = turning_values(a, b);
    vector<long double> zb = turning_values(c, d);
    return path_variation(za, zb);
}

vector<int> primes_in_range(int m, int n) {
    if (m > n) swap(m, n);
    n = max(n, 2);
    vector<bool> is_prime(n + 1, true);
    is_prime[0] = false;
    is_prime[1] = false;
    for (int p = 2; p * p <= n; ++p) {
        if (is_prime[p]) {
            for (int k = p * p; k <= n; k += p) {
                is_prime[k] = false;
            }
        }
    }
    vector<int> primes;
    for (int i = max(2, m); i <= n; ++i) {
        if (is_prime[i] && (i % 2 == 1)) {
            primes.push_back(i);
        }
    }
    return primes;
}

long double compute_G(int m, int n, int threads) {
    vector<int> primes = primes_in_range(m, n);
    vector<pair<int, int>> pairs;
    for (size_t i = 0; i < primes.size(); ++i) {
        for (size_t j = i + 1; j < primes.size(); ++j) {
            pairs.emplace_back(primes[i], primes[j]);
        }
    }

    if (pairs.empty()) {
        return 0.0L;
    }

    if (threads < 1) threads = 1;
    threads = min<int>(threads, static_cast<int>(pairs.size()));

    vector<long double> partial(threads, 0.0L);
    atomic<bool> failed(false);
    if (threads == 1) {
        long double total = 0.0L;
        for (const auto& pq : pairs) {
            int p = pq.first;
            int q = pq.second;
            long double value = compute_F(p, q, p, 2 * q - p);
            if (!isfinite(static_cast<double>(value))) {
                cerr << "Computation failed for pair (" << p << "," << q << ")\n";
                return numeric_limits<long double>::quiet_NaN();
            }
            total += value;
        }
        return total;
    }

    size_t total_pairs = pairs.size();
    size_t chunk = (total_pairs + threads - 1) / threads;
    vector<thread> workers;
    workers.reserve(threads);

    for (int t = 0; t < threads; ++t) {
        size_t start = static_cast<size_t>(t) * chunk;
        size_t end = min(total_pairs, start + chunk);
        workers.emplace_back([&, start, end, t]() {
            long double local = 0.0L;
            for (size_t idx = start; idx < end; ++idx) {
                if (failed.load()) {
                    return;
                }
                int p = pairs[idx].first;
                int q = pairs[idx].second;
                long double value = compute_F(p, q, p, 2 * q - p);
                if (!isfinite(static_cast<double>(value))) {
                    cerr << "Computation failed for pair (" << p << "," << q << ")\n";
                    failed.store(true);
                    return;
                }
                local += value;
            }
            partial[t] = local;
        });
    }

    for (auto& th : workers) {
        th.join();
    }

    long double total = 0.0L;
    for (long double v : partial) {
        total += v;
    }
    if (failed.load()) {
        return numeric_limits<long double>::quiet_NaN();
    }
    return total;
}

bool nearly_equal(long double a, long double b, long double tol) {
    return fabsl(a - b) <= tol;
}

bool run_validation() {
    bool ok = true;
    long double h0 = H(3, 5, 0.0L);
    long double h1 = H(3, 5, 1.0L);
    if (!nearly_equal(h0, 0.0L, 1e-12L) || !nearly_equal(h1, 1.0L, 1e-12L)) {
        cerr << "Validation failed: H(3,5,0/1) expected 0/1, got " << h0 << " / " << h1 << "\n";
        ok = false;
    }

    long double f1 = compute_F(3, 5, 3, 7);
    long double f2 = compute_F(7, 17, 9, 19);
    if (!nearly_equal(f1, 7.01772L, 1e-5L)) {
        cerr << "Validation failed: F(3,5,3,7) " << f1 << "\n";
        ok = false;
    }
    if (!nearly_equal(f2, 26.79578L, 1e-5L)) {
        cerr << "Validation failed: F(7,17,9,19) " << f2 << "\n";
        ok = false;
    }

    long double g = compute_G(3, 20, 1);
    if (!isfinite(static_cast<double>(g))) {
        cerr << "Validation failed: G(3,20) not finite\n";
        ok = false;
    } else if (!nearly_equal(g, 463.80866L, 1e-5L)) {
        cerr << "Validation failed: G(3,20) " << g << "\n";
        ok = false;
    }

    if (ok) {
        cerr << "Validation OK\n";
    }
    return ok;
}

} // namespace

int main(int argc, char** argv) {
    int m = 500;
    int n = 1000;
    int threads = static_cast<int>(thread::hardware_concurrency());
    if (threads == 0) threads = 4;
    bool validate = false;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--validate") {
            validate = true;
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = stoi(argv[++i]);
        } else if (arg == "--range" && i + 2 < argc) {
            m = stoi(argv[++i]);
            n = stoi(argv[++i]);
        } else {
            cerr << "Unknown arg: " << arg << "\n";
            return 1;
        }
    }

    if (validate) {
        if (!run_validation()) {
            return 1;
        }
    }

    long double result = compute_G(m, n, threads);
    cout << fixed << setprecision(5) << static_cast<double>(result) << "\n";
    return 0;
}