#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

// Project Euler 576: Irrational Jumps
//
// For fixed irrational jump length l, the positions are x_k = frac(k*l).
// A gap of length g starting at d catches the walker at the first k with x_k in [d, d+g).
// Equivalently, for fixed x_k, it catches gaps with d in (x_k - g, x_k], intersected with (0, 1-g).
//
// For each l we build the piecewise-constant function T_l(d) = first hit time (number of jumps),
// by "painting" the interval of d-values hit at time k onto the still-uncovered domain, increasing k.
// Then S(l,g,d) = l * T_l(d).
//
// For M(n,g), we need max over d of sum_{p<=n, p prime} S(sqrt(1/p), g, d).
// Each S_p(d) is piecewise-constant; their sum is too. The maximum occurs on one of the induced
// subintervals, so we sweep all segment-start events across primes.

using u64 = std::uint64_t;

static std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    std::vector<int> ps;
    for (int i = 2; i <= n; ++i)
        if (is_prime[i]) ps.push_back(i);
    return ps;
}

struct Seg {
    long double l = 0;
    long double r = 0;
    int k = 0;  // first hit time
};

static std::vector<Seg> build_first_hit_segments(long double alpha, long double g) {
    // Domain for d: (0, 1-g). We'll work on [0, 1-g] and ignore measure-zero endpoint issues.
    const long double domL = 0.0L;
    const long double domR = 1.0L - g;
    const long double eps = 1e-21L;

    std::map<long double, long double> uncovered;
    uncovered.emplace(domL, domR);
    long double remaining = domR - domL;

    std::vector<Seg> segs;
    segs.reserve((size_t)(1.0L / g) + 10);

    long double x = 0.0L;
    int k = 0;
    while (!uncovered.empty() && remaining > eps) {
        ++k;
        x += alpha;
        x -= floorl(x);  // frac

        long double L = x - g;
        long double U = x;
        if (L < domL) L = domL;
        if (U > domR) U = domR;
        if (L + eps >= U) continue;

        // Iterate uncovered intervals that intersect [L, U].
        auto it = uncovered.upper_bound(L);
        if (it != uncovered.begin()) --it;
        while (it != uncovered.end()) {
            const long double a = it->first;
            const long double b = it->second;
            if (b <= L + eps) {
                ++it;
                continue;
            }
            if (a >= U - eps) break;

            const long double ol = std::max(a, L);
            const long double orr = std::min(b, U);
            if (ol + eps < orr) {
                segs.push_back(Seg{ol, orr, k});
                remaining -= (orr - ol);
            }

            // Remove current interval and reinsert leftovers.
            auto it_erase = it++;
            uncovered.erase(it_erase);
            if (a + eps < L) uncovered.emplace(a, std::min(L, b));
            if (U + eps < b) uncovered.emplace(std::max(U, a), b);
        }
    }

    // Sort and merge adjacent segments with the same k (should be rare).
    std::sort(segs.begin(), segs.end(), [](const Seg& A, const Seg& B) { return A.l < B.l; });
    std::vector<Seg> merged;
    merged.reserve(segs.size());
    for (const auto& s : segs) {
        if (merged.empty()) {
            merged.push_back(s);
            continue;
        }
        auto& last = merged.back();
        if (s.k == last.k && fabsl(s.l - last.r) < 1e-18L) {
            last.r = s.r;
        } else {
            merged.push_back(s);
        }
    }

    // Ensure coverage includes the left boundary by filling any tiny gaps (numerical).
    // The exact maximum is unaffected by eps-sized gaps.
    return merged;
}

struct Event {
    long double pos = 0;
    int idx = 0;
    long double new_val = 0;
};

static long double compute_M(int n, long double g) {
    const auto ps = primes_up_to(n);
    const long double domL = 0.0L;
    const long double domR = 1.0L - g;

    const int m = (int)ps.size();
    std::vector<long double> weight(m);
    for (int i = 0; i < m; ++i) weight[i] = sqrtl(1.0L / (long double)ps[i]);

    std::vector<long double> cur(m, 0.0L);
    std::vector<Event> events;
    events.reserve(2000000);

    for (int i = 0; i < m; ++i) {
        const long double alpha = weight[i];
        auto segs = build_first_hit_segments(alpha, g);
        if (segs.empty() || segs.front().l > domL + 1e-15L) {
            std::cerr << "Segment construction failed to cover domain for p=" << ps[i] << "\n";
            std::exit(1);
        }

        // Find the segment covering domL (should be first).
        int j0 = 0;
        while (j0 + 1 < (int)segs.size() && segs[j0].r <= domL) ++j0;
        cur[i] = weight[i] * (long double)segs[j0].k;

        // Create events at subsequent segment starts.
        for (int j = j0 + 1; j < (int)segs.size(); ++j) {
            events.push_back(Event{segs[j].l, i, weight[i] * (long double)segs[j].k});
        }
    }

    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.pos < b.pos; });

    long double sum = 0.0L;
    for (long double v : cur) sum += v;

    long double best = sum;
    long double pos_prev = domL;

    size_t idx = 0;
    while (idx < events.size()) {
        const long double pos = events[idx].pos;
        if (pos > pos_prev + 1e-18L) best = std::max(best, sum);

        // apply all events at this pos
        while (idx < events.size() && fabsl(events[idx].pos - pos) < 1e-18L) {
            const int i = events[idx].idx;
            sum += events[idx].new_val - cur[i];
            cur[i] = events[idx].new_val;
            ++idx;
        }
        pos_prev = pos;
        if (pos_prev > domR) break;
    }
    best = std::max(best, sum);
    return best;
}

int main() {
    {
        const long double got = compute_M(3, 0.06L);
        const long double expected = 29.5425L;
        if (fabsl(got - expected) > 5e-4L) {
            std::cerr << "Validation failed: M(3, 0.06)\n";
            std::cerr << std::fixed << std::setprecision(6) << (double)got << " vs " << (double)expected << "\n";
            return 1;
        }
    }
    {
        const long double got = compute_M(10, 0.01L);
        const long double expected = 266.9010L;
        if (fabsl(got - expected) > 5e-4L) {
            std::cerr << "Validation failed: M(10, 0.01)\n";
            std::cerr << std::fixed << std::setprecision(6) << (double)got << " vs " << (double)expected << "\n";
            return 1;
        }
    }

    const long double ans = compute_M(100, 0.00002L);
    std::cout << std::fixed << std::setprecision(4) << (double)ans << "\n";
    return 0;
}

