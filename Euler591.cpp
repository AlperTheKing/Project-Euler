#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using std::cerr;
using std::cout;
using std::string;
using boost::multiprecision::cpp_dec_float_100;

using Big = cpp_dec_float_100;

struct Approx {
    long long a = 0;
    long long b = 0;
    Big err = 0;
    bool valid = false;
};

static inline long long llabs_ll(long long x) { return x < 0 ? -x : x; }

static long long floor_div(long long a, long long b) {
    __int128 A = static_cast<__int128>(a);
    __int128 B = static_cast<__int128>(b);
    __int128 q = A / B;
    __int128 r = A % B;
    if (r != 0 && A < 0) q -= 1;
    return static_cast<long long>(q);
}

static long long ceil_div(long long a, long long b) {
    __int128 A = static_cast<__int128>(a);
    __int128 B = static_cast<__int128>(b);
    __int128 q = A / B;
    __int128 r = A % B;
    if (r != 0 && A > 0) q += 1;
    return static_cast<long long>(q);
}

static long long floor_ll(const Big& x) {
    long long t = x.convert_to<long long>();
    if (x < 0 && Big(t) != x) return t - 1;
    return t;
}

static long long ceil_ll(const Big& x) {
    long long t = x.convert_to<long long>();
    if (x > 0 && Big(t) != x) return t + 1;
    return t;
}

static long long round_ll(const Big& x) {
    static const Big half("0.5");
    if (x >= 0) return floor_ll(x + half);
    return ceil_ll(x - half);
}

static long long modinv(long long a, long long mod) {
    long long b = mod;
    __int128 x0 = 1, x1 = 0;
    __int128 y0 = 0, y1 = 1;
    long long aa = a;
    while (b != 0) {
        long long q = aa / b;
        long long r = aa % b;
        aa = b;
        b = r;
        __int128 nx = x0 - static_cast<__int128>(q) * x1;
        x0 = x1;
        x1 = nx;
        __int128 ny = y0 - static_cast<__int128>(q) * y1;
        y0 = y1;
        y1 = ny;
    }
    if (aa != 1 && aa != -1) return -1;
    __int128 inv = x0;
    long long res = static_cast<long long>(inv % mod);
    if (res < 0) res += mod;
    return res;
}

static bool is_square(int n) {
    int r = static_cast<int>(std::sqrt(static_cast<double>(n)));
    return r * r == n;
}

static std::vector<int> cf_period_sqrt(int d) {
    int a0 = static_cast<int>(std::sqrt(static_cast<double>(d)));
    if (a0 * a0 == d) return {};
    long long m = 0;
    long long denom = 1;
    long long a = a0;
    std::vector<int> period;
    while (true) {
        m = denom * a - m;
        denom = (d - m * m) / denom;
        a = (a0 + m) / denom;
        period.push_back(static_cast<int>(a));
        if (denom == 1 && m == a0) break;
    }
    return period;
}

static std::vector<std::pair<long long, long long>> convergents_sqrt(int d, long long q_limit) {
    int a0 = static_cast<int>(std::sqrt(static_cast<double>(d)));
    std::vector<int> period = cf_period_sqrt(d);
    std::vector<std::pair<long long, long long>> conv;

    long long p_m2 = 1, q_m2 = 0;
    long long p_m1 = a0, q_m1 = 1;
    conv.push_back({p_m1, q_m1});

    if (period.empty()) return conv;

    size_t i = 0;
    while (true) {
        long long a = period[i % period.size()];
        ++i;
        __int128 p = static_cast<__int128>(a) * p_m1 + p_m2;
        __int128 q = static_cast<__int128>(a) * q_m1 + q_m2;
        long long pn = static_cast<long long>(p);
        long long qn = static_cast<long long>(q);
        conv.push_back({pn, qn});
        if (qn > q_limit) break;
        p_m2 = p_m1;
        q_m2 = q_m1;
        p_m1 = pn;
        q_m1 = qn;
    }
    return conv;
}

static Approx best_BQA_for_d(int d, long long N, int window_j = 25, int window_t = 25) {
    const Big PI = boost::math::constants::pi<Big>();
    const Big TAU = PI - Big(3);

    Big alpha = sqrt(Big(d));

    long long B = floor_ll((Big(N) + PI + Big(2)) / alpha);
    if (B > N) B = N;

    auto conv = convergents_sqrt(d, B);

    Approx best;
    best.err = Big("1e1000");

    auto eval_b = [&](long long b) {
        if (llabs_ll(b) > B || llabs_ll(b) > N) return;
        Big val = alpha * Big(b);
        Big a_real = PI - val;
        long long a0 = round_ll(a_real);

        for (long long da = -1; da <= 1; ++da) {
            long long a = a0 + da;
            if (llabs_ll(a) > N) continue;
            Big err = abs(a_real - Big(a));
            if (!best.valid || err < best.err ||
                (err == best.err && llabs_ll(a) < llabs_ll(best.a)) ||
                (err == best.err && llabs_ll(a) == llabs_ll(best.a) && llabs_ll(b) < llabs_ll(best.b))) {
                best.valid = true;
                best.a = a;
                best.b = b;
                best.err = err;
            }
        }
    };

    for (auto [p_raw, q_raw] : conv) {
        long long q = q_raw;
        long long p = p_raw;
        if (q <= 0) continue;

        long long inv = modinv((p % q + q) % q, q);
        if (inv < 0) continue;

        Big qB = Big(q);
        long long r = (TAU * qB).convert_to<long long>();

        long long j_lo = std::max(0LL, r - static_cast<long long>(window_j));
        long long j_hi = std::min(q - 1, r + static_cast<long long>(window_j));

        Big delta = qB * alpha - Big(p);
        if (delta == 0) continue;

        for (long long j = j_lo; j <= j_hi; ++j) {
            long long b0 = static_cast<long long>((static_cast<__int128>(j) * inv) % q);

            long long t_min = ceil_div(-B - b0, q);
            long long t_max = floor_div(B - b0, q);
            if (t_min > t_max) continue;

            Big diff = Big(j) / qB + (Big(b0) * delta) / qB - TAU;

            Big ratio = (-diff) / delta;
            long long t0 = round_ll(ratio);

            long long tt_lo = std::max(t_min, t0 - static_cast<long long>(window_t));
            long long tt_hi = std::min(t_max, t0 + static_cast<long long>(window_t));
            for (long long t = tt_lo; t <= tt_hi; ++t) {
                long long b = b0 + t * q;
                eval_b(b);
            }

            eval_b(b0 + t_min * q);
            eval_b(b0 + t_max * q);
        }
    }

    return best;
}

static void validation_checkpoints() {
    const long long N13 = 10000000000000LL;

    {
        auto r = best_BQA_for_d(2, 10);
        if (!(r.valid && r.a == 6 && r.b == -2)) {
            cerr << "Validation failed: BQA_2(pi,10) expected 6-2*sqrt(2), got a=" << r.a
                 << " b=" << r.b << "\n";
            std::exit(1);
        }
    }
    {
        auto r = best_BQA_for_d(5, 100);
        if (!(r.valid && r.a == -55 && r.b == 26)) {
            cerr << "Validation failed: BQA_5(pi,100) expected 26*sqrt(5)-55, got a=" << r.a
                 << " b=" << r.b << "\n";
            std::exit(1);
        }
    }
    {
        auto r = best_BQA_for_d(7, 1000000);
        if (!(r.valid && r.a == 560323 && r.b == -211781)) {
            cerr << "Validation failed: BQA_7(pi,1e6) expected 560323-211781*sqrt(7), got a="
                 << r.a << " b=" << r.b << "\n";
            std::exit(1);
        }
    }
    {
        auto r = best_BQA_for_d(2, N13);
        if (!(r.valid && r.a == -6188084046055LL)) {
            cerr << "Validation failed: I_2(BQA_2(pi,1e13)) expected -6188084046055, got a="
                 << r.a << "\n";
            std::exit(1);
        }
    }

    cerr << "All validation checkpoints passed.\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    validation_checkpoints();

    const long long N = 10000000000000LL;

    std::vector<int> ds;
    for (int d = 2; d < 100; ++d) {
        if (is_square(d)) continue;
        ds.push_back(d);
    }

    std::vector<long long> absAs(ds.size(), 0);

    unsigned numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 4;

    std::atomic<size_t> idx{0};
    std::vector<std::thread> threads;
    threads.reserve(numThreads);

    for (unsigned t = 0; t < numThreads; ++t) {
        threads.emplace_back([&]() {
            while (true) {
                size_t i = idx.fetch_add(1, std::memory_order_relaxed);
                if (i >= ds.size()) break;
                int d = ds[i];
                auto r = best_BQA_for_d(d, N);
                if (!r.valid) {
                    cerr << "No valid approximation found for d=" << d << "\n";
                    std::exit(1);
                }
                absAs[i] = llabs_ll(r.a);
            }
        });
    }

    for (auto& th : threads) th.join();

    __int128 sum = 0;
    for (auto v : absAs) sum += v;

    auto print_int128 = [](__int128 x) {
        if (x == 0) {
            cout << '0';
            return;
        }
        if (x < 0) {
            cout << '-';
            x = -x;
        }
        string s;
        while (x > 0) {
            int digit = static_cast<int>(x % 10);
            s.push_back(static_cast<char>('0' + digit));
            x /= 10;
        }
        std::reverse(s.begin(), s.end());
        cout << s;
    };

    print_int128(sum);
    cout << "\n";
    return 0;
}
