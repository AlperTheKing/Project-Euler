#include <cstdint>
#include <iostream>
#include <vector>

using namespace std;

namespace {

constexpr int64_t kMod = 1'000'000'007;

int64_t mod_pow(int64_t base, int64_t exp) {
    int64_t res = 1 % kMod;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1) res = (res * base) % kMod;
        base = (base * base) % kMod;
        exp >>= 1;
    }
    return res;
}

// Compute g(n) using Burnside and O(n) recurrences derived from EGFs.
int64_t compute_g(int n) {
    if (n <= 0) return 0;
    int inv_size = max(n + 2, 9);
    vector<int> inv(inv_size, 0);
    inv[1] = 1;
    for (int i = 2; i < inv_size; ++i) {
        inv[i] = static_cast<int>(
            kMod - (kMod / i) * inv[kMod % i] % kMod);
    }
    int64_t inv2 = inv[2];
    int64_t inv3 = inv[3];
    int64_t inv8 = inv[8];

    int m = n / 2;

    int64_t fact = 1;
    int64_t fact_m = 1;
    for (int i = 1; i <= n; ++i) {
        fact = (fact * i) % kMod;
        if (i == m) fact_m = fact;
    }
    int64_t fact_n = fact;

    // f(n): a_{n+1} = (n a_n + 1/2 a_{n-1}) / (n+1), f(n) = (n!)^2 a_n
    int64_t a_prev = 1; // a_0
    int64_t a_curr = 0; // a_1
    for (int i = 1; i < n; ++i) {
        int64_t a_next = (static_cast<int64_t>(i) * a_curr + inv2 * a_prev) % kMod;
        a_next = (a_next * inv[i + 1]) % kMod;
        a_prev = a_curr;
        a_curr = a_next;
    }
    int64_t a_n = (n == 0) ? a_prev : a_curr;
    int64_t f = fact_n * fact_n % kMod * a_n % kMod;

    // D(n): diagonal reflection counts, b_{n+1} recurrence from EGF.
    int64_t b_n = 0;
    if (n == 0) {
        b_n = 1;
    } else if (n == 1) {
        b_n = 0;
    } else if (n == 2) {
        b_n = inv2;
    } else if (n == 3) {
        b_n = (2 * inv3) % kMod;
    } else {
        int64_t b_nm3 = 1;            // b_0
        int64_t b_nm2 = 0;            // b_1
        int64_t b_nm1 = inv2;         // b_2
        int64_t b_curr = (2 * inv3) % kMod; // b_3
        for (int i = 3; i < n; ++i) {
            int64_t b_next = (2LL * i % kMod * b_curr % kMod
                              - static_cast<int64_t>(i - 2) * b_nm1 % kMod
                              - inv2 * b_nm3 % kMod) % kMod;
            if (b_next < 0) b_next += kMod;
            b_next = (b_next * inv[i + 1]) % kMod;
            b_nm3 = b_nm2;
            b_nm2 = b_nm1;
            b_nm1 = b_curr;
            b_curr = b_next;
        }
        b_n = b_curr;
    }
    int64_t D = fact_n * b_n % kMod;

    // R2(n): rotation 180.
    int64_t R2 = 0;
    if (m > 0) {
        int64_t h_prev = 1; // h_0
        int64_t h_curr = 1; // h_1
        int64_t p_curr = 0; // p_0
        for (int k = 1; k <= m; ++k) {
            p_curr = (4 * p_curr + 2 * h_prev) % kMod; // p_k
            int64_t h_next = ((4LL * k + 1) % kMod * h_curr + 4 * h_prev) % kMod;
            h_next = (h_next * inv[k + 1]) % kMod;
            h_prev = h_curr;
            h_curr = h_next;
        }
        int64_t h_m = h_prev;
        int64_t fact_m_sq = fact_m * fact_m % kMod;
        if (n % 2 == 0) {
            R2 = fact_m_sq * h_m % kMod;
        } else {
            R2 = fact_m_sq * p_curr % kMod;
        }
    }

    // R90(n): rotation 90 (only for even n).
    int64_t R90 = 0;
    if (n % 2 == 0 && m > 0) {
        int64_t c_nm2 = 0;
        int64_t c_nm1 = 0;
        int64_t c_curr = 1; // c_0
        for (int k = 0; k < m; ++k) {
            int64_t c_next = ((2LL * k + 1) % kMod * c_curr - c_nm1 + 2 * c_nm2) % kMod;
            if (c_next < 0) c_next += kMod;
            c_next = (c_next * inv[k + 1]) % kMod;
            c_nm2 = c_nm1;
            c_nm1 = c_curr;
            c_curr = c_next;
        }
        int64_t c_m = c_curr;
        R90 = fact_m * c_m % kMod;
    }

    // V(n): vertical/horizontal reflections.
    int64_t V = 0;
    if (n % 2 == 0) {
        V = fact_n * mod_pow(inv2, n / 2) % kMod;
    }

    int64_t g = (f + R2 + 2 * R90 + 2 * V + 2 * D) % kMod;
    g = g * inv8 % kMod;
    return g;
}

} // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (compute_g(4) != 20) {
        cerr << "Validation failed: g(4) != 20\n";
        return 1;
    }
    if (compute_g(7) != 390816) {
        cerr << "Validation failed: g(7) != 390816\n";
        return 1;
    }
    if (compute_g(8) != 23462347) {
        cerr << "Validation failed: g(8) != 23462347\n";
        return 1;
    }

    int n1 = 1;
    for (int i = 0; i < 7; ++i) n1 *= 7;
    int n2 = 1;
    for (int i = 0; i < 8; ++i) n2 *= 8;

    int64_t g1 = compute_g(n1);
    int64_t g2 = compute_g(n2);
    int64_t ans = (g1 + g2) % kMod;
    cout << ans << "\n";
    return 0;
}
