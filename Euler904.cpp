#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

struct Fraction {
    long long p;
    long long q;
};

struct Candidate {
    bool valid = false;
    long double diff = 0.0L;
    unsigned __int128 area = 0;
    uint64_t sum = 0;
};

static long double tan_theta_from_u(long double u) {
    long double u2 = u * u;
    return 3.0L * u * (u2 - 1.0L) / ((u2 + 1.0L) * (u2 + 1.0L));
}

static Fraction farey_prev(const Fraction& lo, const Fraction& hi, long long N) {
    long long k = (N + hi.q) / lo.q;
    return {k * lo.p - hi.p, k * lo.q - hi.q};
}

static Fraction farey_next(const Fraction& lo, const Fraction& hi, long long N) {
    long long k = (N + lo.q) / hi.q;
    return {k * hi.p - lo.p, k * hi.q - lo.q};
}

static void farey_bounds(long double x, long long N, Fraction& lo, Fraction& hi) {
    std::vector<long long> a;
    long double y = x;
    for (int i = 0; i < 64; ++i) {
        long long ai = static_cast<long long>(std::floor(y));
        a.push_back(ai);
        long double frac = y - static_cast<long double>(ai);
        if (std::fabsl(frac) < 1e-18L) {
            break;
        }
        y = 1.0L / frac;
    }

    std::vector<long long> ps;
    std::vector<long long> qs;
    long long p_minus2 = 0;
    long long p_minus1 = 1;
    long long q_minus2 = 1;
    long long q_minus1 = 0;

    for (long long ai : a) {
        long long p = ai * p_minus1 + p_minus2;
        long long q = ai * q_minus1 + q_minus2;
        ps.push_back(p);
        qs.push_back(q);
        if (q > N) {
            break;
        }
        p_minus2 = p_minus1;
        p_minus1 = p;
        q_minus2 = q_minus1;
        q_minus1 = q;
    }

    int idx = 0;
    while (idx < static_cast<int>(qs.size()) && qs[idx] <= N) {
        ++idx;
    }

    if (idx == 0) {
        lo = {0, 1};
        hi = {1, 0};
        return;
    }

    if (idx == static_cast<int>(qs.size())) {
        lo = {ps.back(), qs.back()};
        hi = lo;
        return;
    }

    long long p_prev = ps[idx - 1];
    long long q_prev = qs[idx - 1];
    long long p_prevprev = (idx >= 2) ? ps[idx - 2] : 1;
    long long q_prevprev = (idx >= 2) ? qs[idx - 2] : 0;

    long long t = (N - q_prevprev) / q_prev;
    long long p_cand = p_prevprev + t * p_prev;
    long long q_cand = q_prevprev + t * q_prev;

    if (static_cast<long double>(p_prev) / static_cast<long double>(q_prev) < x) {
        lo = {p_prev, q_prev};
        hi = {p_cand, q_cand};
    } else {
        lo = {p_cand, q_cand};
        hi = {p_prev, q_prev};
    }
}

static bool is_valid_fraction(const Fraction& f, long long L) {
    if (f.q <= 0 || f.p <= 0) {
        return false;
    }
    if (f.p <= f.q) {
        return false;
    }
    if (((f.p + f.q) & 1LL) == 0) {
        return false;
    }
    if (std::gcd(f.p, f.q) != 1) {
        return false;
    }
    unsigned long long m = static_cast<unsigned long long>(f.p);
    unsigned long long n = static_cast<unsigned long long>(f.q);
    unsigned long long m2 = m * m;
    unsigned long long n2 = n * n;
    return m2 + n2 <= static_cast<unsigned long long>(L);
}

static Fraction find_valid_lower(Fraction lo, Fraction hi, long long N, long long L) {
    if (lo.p <= lo.q) {
        return {0, 0};
    }
    for (int steps = 0; steps < 100000; ++steps) {
        if (is_valid_fraction(lo, L)) {
            return lo;
        }
        Fraction prev = farey_prev(lo, hi, N);
        hi = lo;
        lo = prev;
        if (lo.q == 0) {
            break;
        }
    }
    return {0, 0};
}

static Fraction find_valid_upper(Fraction lo, Fraction hi, long long N, long long L) {
    for (int steps = 0; steps < 100000; ++steps) {
        if (is_valid_fraction(hi, L)) {
            return hi;
        }
        Fraction next = farey_next(lo, hi, N);
        lo = hi;
        hi = next;
        if (hi.q == 0) {
            break;
        }
    }
    return {0, 0};
}

static Candidate evaluate_fraction(const Fraction& f, long double t_alpha, long long L) {
    Candidate cand;
    if (f.p == 0 || f.q == 0) {
        return cand;
    }
    unsigned long long m = static_cast<unsigned long long>(f.p);
    unsigned long long n = static_cast<unsigned long long>(f.q);
    unsigned long long m2 = m * m;
    unsigned long long n2 = n * n;
    unsigned long long a = m2 - n2;
    unsigned long long b = 2ULL * m * n;
    unsigned long long c = m2 + n2;
    unsigned long long k = static_cast<unsigned long long>(L) / c;
    if (k == 0) {
        return cand;
    }

    long double t = 3.0L * static_cast<long double>(a) * static_cast<long double>(b)
        / (2.0L * static_cast<long double>(c) * static_cast<long double>(c));
    long double diff = std::fabsl(t - t_alpha) / (1.0L + t * t_alpha);

    cand.valid = true;
    cand.diff = diff;
    cand.sum = static_cast<uint64_t>(k) * static_cast<uint64_t>(a + b + c);
    cand.area = static_cast<unsigned __int128>(k) * static_cast<unsigned __int128>(k)
        * static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b);
    return cand;
}

static Candidate best_candidate_for_u(long double u_target, long double t_alpha, long long L) {
    Candidate best;
    if (u_target <= 1.0L) {
        return best;
    }

    long double R = sqrtl(static_cast<long double>(L));
    long double denom = sqrtl(1.0L + u_target * u_target);
    long long N = static_cast<long long>(std::floor(R / denom));
    if (N < 1) {
        N = 1;
    }

    Fraction lo{0, 1};
    Fraction hi{1, 0};
    farey_bounds(u_target, N, lo, hi);

    Fraction lo_valid = find_valid_lower(lo, hi, N, L);
    Fraction hi_valid = find_valid_upper(lo, hi, N, L);

    if (lo_valid.q != 0) {
        best = evaluate_fraction(lo_valid, t_alpha, L);
    }
    if (hi_valid.q != 0) {
        Candidate cand = evaluate_fraction(hi_valid, t_alpha, L);
        if (!best.valid) {
            best = cand;
        } else if (cand.valid) {
            const long double eps = 1e-21L;
            if (cand.diff + eps < best.diff) {
                best = cand;
            } else if (std::fabsl(cand.diff - best.diff) <= eps) {
                if (cand.area > best.area) {
                    best = cand;
                }
            }
        }
    }

    return best;
}

static uint64_t compute_f(long double alpha_deg, long long L) {
    static const long double kPi = acosl(-1.0L);
    static const long double kUMin = 1.0L + sqrtl(2.0L);
    static const long double kTMax = 0.75L;

    long double alpha_rad = alpha_deg * kPi / 180.0L;
    long double t_alpha = tanl(alpha_rad);

    long double u_low = 1.0L;
    long double u_high = kUMin;
    if (t_alpha < kTMax) {
        long double low = 1.0L;
        long double high = kUMin;
        for (int iter = 0; iter < 100; ++iter) {
            long double mid = 0.5L * (low + high);
            if (tan_theta_from_u(mid) < t_alpha) {
                low = mid;
            } else {
                high = mid;
            }
        }
        u_low = 0.5L * (low + high);

        low = kUMin;
        high = std::max(kUMin + 1.0L, 3.0L / t_alpha + 2.0L);
        while (tan_theta_from_u(high) > t_alpha) {
            high *= 2.0L;
        }
        for (int iter = 0; iter < 100; ++iter) {
            long double mid = 0.5L * (low + high);
            if (tan_theta_from_u(mid) > t_alpha) {
                low = mid;
            } else {
                high = mid;
            }
        }
        u_high = 0.5L * (low + high);
    }

    Candidate best = best_candidate_for_u(u_low, t_alpha, L);
    Candidate alt = best_candidate_for_u(u_high, t_alpha, L);
    if (!best.valid || (alt.valid && (alt.diff < best.diff - 1e-21L
        || (std::fabsl(alt.diff - best.diff) <= 1e-21L && alt.area > best.area)))) {
        best = alt;
    }

    return best.sum;
}

static uint64_t compute_F(int N, long long L) {
    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 1;
    }

    std::vector<uint64_t> partial(threads, 0);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned int t = 0; t < threads; ++t) {
        pool.emplace_back([t, threads, N, L, &partial]() {
            uint64_t local = 0;
            for (int i = 1 + static_cast<int>(t); i <= N; i += static_cast<int>(threads)) {
                long double alpha = cbrtl(static_cast<long double>(i));
                local += compute_f(alpha, L);
            }
            partial[t] = local;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    uint64_t total = 0;
    for (uint64_t v : partial) {
        total += v;
    }
    return total;
}

static bool run_validation() {
    bool ok = true;
    if (compute_f(30.0L, 100LL) != 198ULL) {
        std::cerr << "Validation failed: f(30, 10^2) != 198\n";
        ok = false;
    }
    if (compute_f(10.0L, 1000000LL) != 1600158ULL) {
        std::cerr << "Validation failed: f(10, 10^6) != 1600158\n";
        ok = false;
    }
    if (compute_F(10, 1000000LL) != 16684370ULL) {
        std::cerr << "Validation failed: F(10, 10^6) != 16684370\n";
        ok = false;
    }
    if (ok) {
        std::cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main() {
    if (!run_validation()) {
        return 1;
    }

    uint64_t result = compute_F(45000, 10000000000LL);
    std::cout << result << "\n";
    return 0;
}
