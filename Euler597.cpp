#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

using BigFloat = boost::multiprecision::cpp_dec_float_50;

namespace {

struct Table {
    int max_n;
    int max_m;
    vector<vector<BigFloat>> S;

    Table(int n, int m) : max_n(n), max_m(m), S(n + 1, vector<BigFloat>(m + 1, BigFloat(0))) {
        build();
    }

    void build() {
        for (int m = 0; m <= max_m; ++m) {
            S[0][m] = BigFloat(1);
            S[1][m] = BigFloat(1);
        }

        for (int n = 2; n <= max_n; ++n) {
            int m0 = n - 1;
            if (m0 <= max_m) {
                BigFloat sign = ((n - 1) & 1) ? BigFloat(-1) : BigFloat(1);
                S[n][m0] = sign * S[n - 1][m0];
            }

            for (int m = n; m <= max_m; ++m) {
                long long D_int = 1LL * n * m - 1LL * n * (n - 1) / 2;
                BigFloat D = BigFloat(D_int);
                BigFloat total = 0;
                for (int k = 1; k <= n; ++k) {
                    int d_k = m - (k - 1);
                    BigFloat weight = BigFloat(d_k) / D;
                    total += weight * S[k][k - 1] * S[n - k][m - k];
                }
                S[n][m] = total;
            }
        }
    }

    BigFloat parity_sign(int n, int m) const {
        return S[n][m];
    }

    BigFloat even_prob(int n, int m) const {
        return (S[n][m] + BigFloat(1)) / BigFloat(2);
    }
};

bool approx_equal(const BigFloat& a, const BigFloat& b, const BigFloat& eps) {
    using boost::multiprecision::abs;
    return abs(a - b) <= eps;
}

bool validate(const Table& table) {
    const BigFloat eps = BigFloat("1e-30");
    bool ok = true;

    const BigFloat p3 = table.even_prob(3, 4);
    const BigFloat p3_expected = BigFloat(56) / BigFloat(135);
    if (!approx_equal(p3, p3_expected, eps)) {
        cerr << "Validation failed: p(3,160) got " << p3
             << " expected " << p3_expected << "\n";
        ok = false;
    }

    const BigFloat p4 = table.even_prob(4, 10);
    const BigFloat p4_expected = BigFloat(521) / BigFloat(1020);
    if (!approx_equal(p4, p4_expected, eps)) {
        cerr << "Validation failed: p(4,400) got " << p4
             << " expected " << p4_expected << "\n";
        ok = false;
    }

    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

}  // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int n = 13;
    const int L = 1800;
    const int gap = 40;
    if (L % gap != 0) {
        cerr << "Expected L to be a multiple of 40.\n";
        return 1;
    }
    const int m = L / gap;

    Table table(n, m);
    if (!validate(table)) {
        return 1;
    }

    const BigFloat ans = table.even_prob(n, m);
    if (ans < BigFloat(0) || ans > BigFloat(1)) {
        cerr << "Probability out of range: " << ans << "\n";
        return 1;
    }

    cout << fixed << setprecision(10) << ans << "\n";
    return 0;
}
