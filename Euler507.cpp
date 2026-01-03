#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr int64_t MOD = 10000000LL;

struct Vec3 {
    int64_t x = 0;
    int64_t y = 0;
    int64_t z = 0;
};

inline __int128 abs128(int64_t v) {
    return v < 0 ? -static_cast<__int128>(v) : static_cast<__int128>(v);
}

inline __int128 l1(const Vec3& v) {
    return abs128(v.x) + abs128(v.y) + abs128(v.z);
}

int64_t floor_div(int64_t a, int64_t b) {
    int64_t q = a / b;
    int64_t r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) --q;
    return q;
}

__int128 l1_diff(const Vec3& w, const Vec3& v, int64_t t) {
    __int128 sum = 0;
    __int128 val = static_cast<__int128>(w.x) - static_cast<__int128>(t) * v.x;
    sum += val < 0 ? -val : val;
    val = static_cast<__int128>(w.y) - static_cast<__int128>(t) * v.y;
    sum += val < 0 ? -val : val;
    val = static_cast<__int128>(w.z) - static_cast<__int128>(t) * v.z;
    sum += val < 0 ? -val : val;
    return sum;
}

int64_t best_t(const Vec3& v, const Vec3& w) {
    int64_t cand[8];
    int count = 0;
    cand[count++] = 0;

    auto add = [&](int64_t t) { cand[count++] = t; };
    if (v.x != 0) {
        int64_t q = floor_div(w.x, v.x);
        add(q);
        add(q + 1);
    }
    if (v.y != 0) {
        int64_t q = floor_div(w.y, v.y);
        add(q);
        add(q + 1);
    }
    if (v.z != 0) {
        int64_t q = floor_div(w.z, v.z);
        add(q);
        add(q + 1);
    }

    int64_t best = cand[0];
    __int128 best_val = l1_diff(w, v, best);
    for (int i = 1; i < count; ++i) {
        __int128 val = l1_diff(w, v, cand[i]);
        if (val < best_val) {
            best_val = val;
            best = cand[i];
        }
    }
    return best;
}

__int128 shortest_l1(Vec3 a, Vec3 b) {
    while (true) {
        if (l1(b) < l1(a)) swap(a, b);
        int64_t t = best_t(a, b);
        Vec3 nb{
            static_cast<int64_t>(static_cast<__int128>(b.x) - static_cast<__int128>(t) * a.x),
            static_cast<int64_t>(static_cast<__int128>(b.y) - static_cast<__int128>(t) * a.y),
            static_cast<int64_t>(static_cast<__int128>(b.z) - static_cast<__int128>(t) * a.z),
        };
        if (l1(nb) < l1(b)) {
            b = nb;
        } else {
            break;
        }
    }

    __int128 la = l1(a);
    __int128 lb = l1(b);
    if (la == 0) return lb;
    if (lb == 0) return la;
    return la < lb ? la : lb;
}

struct Mat3 {
    int64_t m[3][3];
};

Mat3 mul(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            __int128 sum = 0;
            for (int k = 0; k < 3; ++k) {
                sum += static_cast<__int128>(a.m[i][k]) * b.m[k][j];
            }
            r.m[i][j] = static_cast<int64_t>(sum % MOD);
        }
    }
    return r;
}

Mat3 mat_pow(Mat3 base, int64_t exp) {
    Mat3 res{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) res.m[i][j] = (i == j) ? 1 : 0;
    }
    while (exp > 0) {
        if (exp & 1LL) res = mul(res, base);
        base = mul(base, base);
        exp >>= 1;
    }
    return res;
}

struct Triple {
    int64_t tn = 0;
    int64_t tn1 = 0;
    int64_t tn2 = 0;
};

Triple trib_state(int64_t n) {
    if (n == 1) return {0, 0, 1};
    if (n == 2) return {1, 0, 0};
    Mat3 base{{
        {1, 1, 1},
        {1, 0, 0},
        {0, 1, 0},
    }};
    Mat3 p = mat_pow(base, n - 2);
    return {p.m[0][0], p.m[1][0], p.m[2][0]};
}

__int128 sum_range(int64_t n_start, int64_t n_end) {
    if (n_start > n_end) return 0;
    int64_t idx_start = 12 * n_start - 11;
    int64_t idx_end = 12 * n_end;

    Triple state = trib_state(idx_start);
    int64_t a = state.tn2;
    int64_t b = state.tn1;
    int64_t c = state.tn;

    int64_t total = idx_end - idx_start + 1;
    int64_t rvals[12];
    int pos = 0;
    __int128 sum = 0;

    for (int64_t i = 0; i < total; ++i) {
        rvals[pos++] = c;
        if (pos == 12) {
            Vec3 v{
                rvals[0] - rvals[1],
                rvals[2] + rvals[3],
                rvals[4] * rvals[5],
            };
            Vec3 w{
                rvals[6] - rvals[7],
                rvals[8] + rvals[9],
                rvals[10] * rvals[11],
            };
            sum += shortest_l1(v, w);
            pos = 0;
        }

        int64_t next = a + b + c;
        if (next >= MOD) {
            next -= MOD;
            if (next >= MOD) next -= MOD;
        }
        a = b;
        b = c;
        c = next;
    }

    return sum;
}

__int128 compute_sum(int64_t n, unsigned threads) {
    if (n <= 0) return 0;
    if (threads == 0) threads = 1;
    if (n < 1000) threads = 1;
    if (threads > static_cast<unsigned>(n)) threads = static_cast<unsigned>(n);

    vector<__int128> partial(threads, 0);
    vector<thread> pool;
    pool.reserve(threads);

    int64_t base = n / threads;
    int64_t rem = n % threads;
    int64_t start = 1;

    for (unsigned i = 0; i < threads; ++i) {
        int64_t count = base + (i < static_cast<unsigned>(rem) ? 1 : 0);
        int64_t n_start = start;
        int64_t n_end = start + count - 1;
        start = n_end + 1;

        pool.emplace_back([&, i, n_start, n_end, count]() {
            if (count == 0) {
                partial[i] = 0;
                return;
            }
            partial[i] = sum_range(n_start, n_end);
        });
    }

    for (auto& th : pool) th.join();

    __int128 total = 0;
    for (const auto& v : partial) total += v;
    return total;
}

string to_string128(__int128 value) {
    if (value == 0) return "0";
    bool neg = value < 0;
    if (neg) value = -value;
    string out;
    while (value > 0) {
        int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (neg) out.push_back('-');
    reverse(out.begin(), out.end());
    return out;
}

bool validate() {
    int64_t r[13];
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    for (int i = 3; i <= 12; ++i) {
        int64_t next = r[i - 1] + r[i - 2] + r[i - 3];
        if (next >= MOD) {
            next -= MOD;
            if (next >= MOD) next -= MOD;
        }
        r[i] = next;
    }

    Vec3 v{r[1] - r[2], r[3] + r[4], r[5] * r[6]};
    Vec3 w{r[7] - r[8], r[9] + r[10], r[11] * r[12]};
    if (!(v.x == -1 && v.y == 3 && v.z == 28 && w.x == -11 && w.y == 125 && w.z == 40826)) {
        cerr << "Validation failed: first vector pair mismatch." << '\n';
        return false;
    }

    __int128 s1 = compute_sum(1, 1);
    if (s1 != 32) {
        cerr << "Validation failed: S(1)=" << to_string128(s1) << ", expected 32." << '\n';
        return false;
    }

    __int128 s10 = compute_sum(10, 1);
    if (s10 != 130762273722LL) {
        cerr << "Validation failed: sum S(1..10)=" << to_string128(s10)
             << ", expected 130762273722." << '\n';
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (!validate()) return 1;

    int64_t n = 20000000;
    if (argc > 1) n = max<int64_t>(1, atoll(argv[1]));
    unsigned threads = thread::hardware_concurrency();
    if (argc > 2) threads = max(1, atoi(argv[2]));

    __int128 result = compute_sum(n, threads);
    cout << to_string128(result) << '\n';
    return 0;
}
