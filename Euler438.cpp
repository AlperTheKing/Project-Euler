#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;
using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

constexpr double kEps = 1e-9;
constexpr int kN = 7;
constexpr int kM = 2;
constexpr int kLower[kN + 1] = {1, 28, 322, 1960, 6769, 13132, 13068, 5040};
constexpr int kUpper[kN + 1] = {1, 35, 511, 4025, 18424, 48860, 69264, 40320};

struct System {
    Mat lower;
    Mat upper;
    Mat equal;

    void add_constraint(Vec c) {
        double coef = -c.back();
        c.pop_back();
        if (std::abs(coef) < kEps) {
            equal.push_back(std::move(c));
            return;
        }
        for (double& v : c) v /= coef;
        if (coef > kEps) {
            lower.push_back(std::move(c));
        } else {
            upper.push_back(std::move(c));
        }
    }
};

Vec diff_vec(const Vec& a, const Vec& b) {
    Vec c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) c[i] = a[i] - b[i];
    return c;
}

System eliminate_one(const System& s) {
    System out;
    for (const Vec& e : s.equal) out.add_constraint(e);
    for (const Vec& lo : s.lower) {
        for (const Vec& hi : s.upper) {
            out.add_constraint(diff_vec(lo, hi));
        }
    }
    return out;
}

double eval_linear(const Vec& a, const std::vector<int>& x) {
    double v = a[0];
    for (std::size_t i = 0; i < x.size(); ++i) v += a[i + 1] * static_cast<double>(x[i]);
    return v;
}

int tight_floor(const double x) {
    return static_cast<int>(std::floor(x + kEps));
}

int tight_ceil(const double x) {
    return static_cast<int>(std::ceil(x - kEps));
}

std::pair<int, int> infer_range(const System& s, const std::vector<int>& x) {
    for (const Vec& e : s.equal) {
        if (eval_linear(e, x) > kEps) return {1, 0};
    }
    int l = -1000000000;
    int r = 1000000000;
    for (const Vec& lo : s.lower) l = std::max(l, tight_ceil(eval_linear(lo, x)));
    for (const Vec& hi : s.upper) r = std::min(r, tight_floor(eval_linear(hi, x)));
    return {l, r};
}

bool valid_endpoint(std::vector<int> x, const int last) {
    x.push_back(last);
    for (int k = 1; k <= kN + 1; ++k) {
        int val = 1;
        for (int i = 0; i < kN; ++i) {
            val = k * val + ((i % 2 == 1) ? x[i] : -x[i]);
        }
        if (val == 0) {
            int der = kN;
            for (int i = 0; i < kN - 1; ++i) {
                der = k * der + (kN - i - 1) * ((i % 2 == 1) ? x[i] : -x[i]);
            }
            if (((kN + k) % 2 == 1) ? (der >= 0) : (der <= 0)) return false;
        }
    }
    return true;
}

std::string to_string_i128(i128 v) {
    if (v == 0) return "0";
    bool neg = (v < 0);
    if (neg) v = -v;
    std::string s;
    while (v > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(v % 10)));
        v /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

struct Solver {
    std::vector<System> systems;
    std::vector<int> x;
    i128 answer = 0;

    Solver() : systems(kN + 2) {
        Mat base;
        base.reserve(kN + 1);
        for (int k = 1; k <= kN + 1; ++k) {
            Vec d(kN + 1, 0.0);
            d[kN] = 1.0;
            for (int j = kN - 1; j >= 0; --j) d[j] = -k * d[j + 1];
            if (k % 2 == 1) {
                for (double& v : d) v = -v;
            }
            base.push_back(std::move(d));
        }
        systems[kN + 1] = System{{}, {}, base};
        for (int k = kN; k > kM; --k) {
            systems[k] = eliminate_one(systems[k + 1]);
        }
    }

    void dfs(const int idx) {
        int l = kLower[idx];
        int r = kUpper[idx] - 1;
        if (idx > kM) {
            const auto p = infer_range(systems[idx], x);
            l = std::max(l, p.first);
            r = std::min(r, p.second);
        }
        if (l > r) return;

        if (idx == kN) {
            while (l <= r && !valid_endpoint(x, l)) ++l;
            while (l <= r && !valid_endpoint(x, r)) --r;
            if (l > r) return;

            const i64 count = static_cast<i64>(r - l + 1);
            const i64 linear_sum = static_cast<i64>(count * static_cast<i64>(l + r) / 2);
            i64 prefix_sum = 0;
            for (int v : x) prefix_sum += v;
            answer += static_cast<i128>(linear_sum) + static_cast<i128>(prefix_sum) * count;
            return;
        }

        for (int v = l; v <= r; ++v) {
            x.push_back(v);
            dfs(idx + 1);
            x.pop_back();
        }
    }

    i128 solve() {
        dfs(1);
        return answer;
    }
};

}  // namespace

int main() {
    Solver solver;
    std::cout << to_string_i128(solver.solve()) << '\n';
    return 0;
}
