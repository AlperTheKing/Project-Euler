#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;

namespace {

using BigFloat = boost::multiprecision::cpp_dec_float_50;
using Poly = vector<BigFloat>;

struct Result {
    BigFloat prob_zero;
    BigFloat total_mass;
};

static Result compute_probability(int n) {
    assert(n > 1);
    const int m = n - 1;
    const int max_deg = 2 * m;
    const int offset = m;

    auto make_zero_poly = [&]() { return Poly(max_deg + 1, BigFloat(0)); };

    vector<Poly> curr(2 * m + 1, make_zero_poly());
    curr[offset][0] = BigFloat(1);

    for (int step = 1; step <= m; ++step) {
        vector<Poly> next(2 * m + 1, make_zero_poly());

        for (int s = -(step - 1); s <= (step - 1); ++s) {
            const Poly& P = curr[offset + s];

            BigFloat total = 0;
            BigFloat total_r = 0;
            for (int i = 0; i <= max_deg; ++i) {
                if (P[i] == 0) continue;
                total += P[i] / BigFloat(i + 1);
                total_r += P[i] / BigFloat(i + 2);
            }

            Poly A(max_deg + 1, BigFloat(0));
            Poly B(max_deg + 1, BigFloat(0));
            for (int i = 0; i <= max_deg; ++i) {
                if (P[i] == 0) continue;
                if (i + 1 <= max_deg) A[i + 1] = P[i] / BigFloat(i + 1);
                if (i + 2 <= max_deg) B[i + 2] = P[i] / BigFloat(i + 2);
            }

            Poly TminusA(max_deg + 1, BigFloat(0));
            TminusA[0] = total;
            for (int i = 1; i <= max_deg; ++i) TminusA[i] = -A[i];

            Poly& same = next[offset + s];
            Poly& plus = next[offset + s + 1];
            Poly& minus = next[offset + s - 1];

            for (int i = 0; i <= max_deg; ++i) {
                BigFloat i1 = A[i] + B[i];
                if (i > 0) i1 -= A[i - 1];

                BigFloat i2 = TminusA[i] + B[i];
                if (i > 0) i2 += TminusA[i - 1];
                if (i == 0) i2 -= total_r;

                same[i] += i1 + i2;

                BigFloat jplus = -B[i];
                if (i > 0) jplus += A[i - 1];
                plus[i] += jplus;

                BigFloat jminus = -B[i];
                if (i == 0) jminus += total_r;
                if (i > 0) jminus -= TminusA[i - 1];
                minus[i] += jminus;
            }
        }

        curr.swap(next);
    }

    auto integrate_0_1 = [&](const Poly& P) {
        BigFloat res = 0;
        for (int i = 0; i <= max_deg; ++i) {
            if (P[i] == 0) continue;
            res += P[i] / BigFloat(i + 1);
        }
        return res;
    };

    Result out;
    out.prob_zero = integrate_0_1(curr[offset]);

    BigFloat total_mass = 0;
    for (int s = -m; s <= m; ++s) {
        total_mass += integrate_0_1(curr[offset + s]);
    }
    out.total_mass = total_mass;
    return out;
}

static void validate() {
    using boost::multiprecision::abs;
    const BigFloat eps = BigFloat("1e-22");

    const Result r3 = compute_probability(3);
    const Result r5 = compute_probability(5);

    const BigFloat expected3 = BigFloat(11) / BigFloat(20);
    const BigFloat expected5 = BigFloat(15619) / BigFloat(36288);

    if (abs(r3.prob_zero - expected3) > eps) {
        cerr << "Validation failed for P(3)." << '\n';
        std::exit(1);
    }
    if (abs(r5.prob_zero - expected5) > eps) {
        cerr << "Validation failed for P(5)." << '\n';
        std::exit(1);
    }
    if (abs(r5.total_mass - BigFloat(1)) > eps) {
        cerr << "Validation failed for normalization." << '\n';
        std::exit(1);
    }
}

}

int main() {
    using boost::multiprecision::abs;
    validate();

    const Result r80 = compute_probability(80);
    const BigFloat eps = BigFloat("1e-22");
    if (abs(r80.total_mass - BigFloat(1)) > eps) {
        cerr << "Normalization check failed for n=80." << '\n';
        return 1;
    }

    cout << fixed << setprecision(10) << r80.prob_zero << '\n';
    return 0;
}
