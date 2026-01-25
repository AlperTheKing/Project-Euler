#include <cstdint>
#include <iostream>
#include <vector>

using namespace std;

static constexpr uint64_t MOD = 888888883ULL;

static uint64_t mod_pow(uint64_t base, uint64_t exp) {
    uint64_t result = 1 % MOD;
    base %= MOD;
    while (exp > 0) {
        if (exp & 1ULL) result = (result * base) % MOD;
        base = (base * base) % MOD;
        exp >>= 1ULL;
    }
    return result;
}

static vector<uint64_t> build_factorials(int max_n, vector<uint64_t>& invfac) {
    vector<uint64_t> fac(max_n + 1, 1);
    for (int i = 1; i <= max_n; ++i) {
        fac[i] = (fac[i - 1] * static_cast<uint64_t>(i)) % MOD;
    }
    invfac.assign(max_n + 1, 1);
    invfac[max_n] = mod_pow(fac[max_n], MOD - 2);
    for (int i = max_n; i >= 1; --i) {
        invfac[i - 1] = (invfac[i] * static_cast<uint64_t>(i)) % MOD;
    }
    return fac;
}

static uint64_t multinomial_mod(int n, int a, int b, int c,
                                const vector<uint64_t>& fac,
                                const vector<uint64_t>& invfac) {
    uint64_t res = fac[n];
    res = (res * invfac[a]) % MOD;
    res = (res * invfac[b]) % MOD;
    res = (res * invfac[c]) % MOD;
    return res;
}

static uint64_t compute_N_mod(int X, int Y, int Z,
                              const vector<uint64_t>& fac,
                              const vector<uint64_t>& invfac) {
    int px = X & 1;
    int py = Y & 1;
    int pz = Z & 1;
    if (px != py || py != pz) return 0;

    int S = X + Y + Z;
    // Neutral strings satisfy inversion parity == ceil(S/2) mod 2.
    // Let total be the multinomial count and diff = even - odd.
    // For all even counts, diff equals the multinomial of half counts.
    // For all odd counts, diff is zero so both parities have total/2 strings.
    uint64_t total = multinomial_mod(S, X, Y, Z, fac, invfac);
    uint64_t inv2 = (MOD + 1) / 2;

    if (px == 1) {
        // All odd counts -> even/odd inversion counts are equal.
        return (total * inv2) % MOD;
    }

    // All even counts -> difference equals multinomial of half counts.
    uint64_t half = multinomial_mod(S / 2, X / 2, Y / 2, Z / 2, fac, invfac);
    if ((S / 2) & 1) {
        uint64_t diff = (total + MOD - half) % MOD;
        return (diff * inv2) % MOD;
    }
    uint64_t sum = total + half;
    if (sum >= MOD) sum -= MOD;
    return (sum * inv2) % MOD;
}

static __int128 factorial128(int n) {
    __int128 res = 1;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

static unsigned long long multinomial128(int n, int a, int b, int c) {
    __int128 res = factorial128(n);
    res /= factorial128(a);
    res /= factorial128(b);
    res /= factorial128(c);
    return static_cast<unsigned long long>(res);
}

static unsigned long long compute_N_exact_small(int X, int Y, int Z) {
    int px = X & 1;
    int py = Y & 1;
    int pz = Z & 1;
    if (px != py || py != pz) return 0;

    int S = X + Y + Z;
    unsigned long long total = multinomial128(S, X, Y, Z);

    if (px == 1) {
        return total / 2;
    }

    unsigned long long half = multinomial128(S / 2, X / 2, Y / 2, Z / 2);
    if ((S / 2) & 1) {
        return (total - half) / 2;
    }
    return (total + half) / 2;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int max_i = 87;
    vector<int> cubes(max_i + 1, 0);
    for (int i = 0; i <= max_i; ++i) {
        cubes[i] = i * i * i;
    }
    int max_sum = 3 * cubes[max_i];

    vector<uint64_t> invfac;
    vector<uint64_t> fac = build_factorials(max_sum, invfac);

    // Validation checkpoints from the statement.
    {
        unsigned long long n222 = compute_N_exact_small(2, 2, 2);
        unsigned long long n888 = compute_N_exact_small(8, 8, 8);
        if (n222 != 42ULL || n888 != 4732773210ULL) {
            cerr << "Validation failed: N(2,2,2)=" << n222
                 << ", N(8,8,8)=" << n888 << "\n";
            return 1;
        }
    }

    uint64_t answer = 0;
    for (int i = 0; i <= max_i; ++i) {
        int X = cubes[i];
        for (int j = 0; j <= max_i; ++j) {
            int Y = cubes[j];
            for (int k = 0; k <= max_i; ++k) {
                int Z = cubes[k];
                answer += compute_N_mod(X, Y, Z, fac, invfac);
                if (answer >= MOD) answer -= MOD;
            }
        }
    }

    cout << answer % MOD << "\n";
    return 0;
}
