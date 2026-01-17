#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

using int64 = long long;

namespace {

constexpr int64 MOD = 999999017LL;
constexpr int64 DEFAULT_N = 99999999019LL;

inline int64 mod_norm(int64 x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return x;
}

inline int64 mod_add(int64 a, int64 b) {
    a += b;
    if (a >= MOD) a -= MOD;
    if (a < 0) a += MOD;
    return a;
}

inline int64 mod_sub(int64 a, int64 b) {
    a -= b;
    if (a < 0) a += MOD;
    return a;
}

inline int64 mod_mul(int64 a, int64 b) {
    return static_cast<int64>((static_cast<__int128>(mod_norm(a)) * mod_norm(b)) % MOD);
}

int64 egcd(int64 a, int64 b, int64& x, int64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int64 x1, y1;
    int64 g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return g;
}

int64 mod_inv(int64 a) {
    int64 x, y;
    int64 g = egcd(a, MOD, x, y);
    if (g != 1) return 0;
    return mod_norm(x);
}

constexpr int64 INV2 = (MOD + 1) / 2;
const int64 INV6 = mod_inv(6);

int64 sum_arith(int64 l, int64 r) {
    if (l > r) return 0;
    int64 cnt = (r - l + 1) % MOD;
    int64 s = mod_norm(l + r);
    return mod_mul(mod_mul(s, cnt), INV2);
}

int64 sum_sq(int64 n) {
    n = mod_norm(n);
    int64 a = n;
    int64 b = mod_norm(n + 1);
    int64 c = mod_norm(2 * n + 1);
    return mod_mul(mod_mul(mod_mul(a, b), c), INV6);
}

struct Solver {
    int64 N;
    int LIM = 0;

    std::vector<int> primes;
    std::vector<int8_t> mu;
    std::vector<uint32_t> phi;
    std::vector<uint8_t> is_comp;
    std::vector<int32_t> prefH;
    std::vector<int32_t> prefG;

    std::unordered_map<int64, int32_t> memoH;
    std::unordered_map<int64, int32_t> memoG;

    explicit Solver(int64 n) : N(n) {
        long double x = std::pow(static_cast<long double>(N), 2.0L / 3.0L);
        LIM = static_cast<int>(x + 10);
        if (LIM < 100) LIM = 100;

        sieve();

        memoH.reserve(1 << 20);
        memoG.reserve(1 << 20);
        memoH.max_load_factor(0.7f);
        memoG.max_load_factor(0.7f);
    }

    void sieve() {
        mu.assign(LIM + 1, 0);
        phi.assign(LIM + 1, 0);
        is_comp.assign(LIM + 1, 0);
        prefH.assign(LIM + 1, 0);
        prefG.assign(LIM + 1, 0);

        primes.clear();
        primes.reserve(LIM / 10);

        mu[1] = 1;
        phi[1] = 1;

        for (int i = 2; i <= LIM; ++i) {
            if (!is_comp[i]) {
                primes.push_back(i);
                mu[i] = -1;
                phi[i] = static_cast<uint32_t>(i - 1);
            }
            for (int p : primes) {
                int64 v = static_cast<int64>(i) * p;
                if (v > LIM) break;
                is_comp[static_cast<size_t>(v)] = 1;
                if (i % p == 0) {
                    mu[static_cast<size_t>(v)] = 0;
                    phi[static_cast<size_t>(v)] = phi[i] * static_cast<uint32_t>(p);
                    break;
                }
                mu[static_cast<size_t>(v)] = static_cast<int8_t>(-mu[i]);
                phi[static_cast<size_t>(v)] = phi[i] * static_cast<uint32_t>(p - 1);
            }
        }

        for (int i = 1; i <= LIM; ++i) {
            int64 addH = mod_mul(static_cast<int64>(mu[i]), static_cast<int64>(i));
            prefH[i] = static_cast<int32_t>(mod_add(prefH[i - 1], addH));

            int64 addG = mod_mul(static_cast<int64>(i), static_cast<int64>(phi[i]));
            prefG[i] = static_cast<int32_t>(mod_add(prefG[i - 1], addG));
        }
    }

    int32_t H(int64 n) {
        if (n <= 0) return 0;
        if (n <= LIM) return prefH[static_cast<size_t>(n)];
        auto it = memoH.find(n);
        if (it != memoH.end()) return it->second;

        int64 res = 1 % MOD;
        int64 l = 2;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            int64 coef = sum_arith(l, r);
            res = mod_sub(res, mod_mul(coef, H(v)));
            l = r + 1;
        }

        int32_t out = static_cast<int32_t>(res);
        memoH.emplace(n, out);
        return out;
    }

    int32_t G(int64 n) {
        if (n <= 0) return 0;
        if (n <= LIM) return prefG[static_cast<size_t>(n)];
        auto it = memoG.find(n);
        if (it != memoG.end()) return it->second;

        int64 res = 0;
        int64 l = 1;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            int64 mu_seg = mod_sub(H(r), H(l - 1));
            res = mod_add(res, mod_mul(mu_seg, sum_sq(v)));
            l = r + 1;
        }

        int32_t out = static_cast<int32_t>(res);
        memoG.emplace(n, out);
        return out;
    }

    int64 T(int64 n) {
        int64 res = 0;
        int64 l = 1;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            int64 seg = mod_sub(G(r), G(l - 1));
            res = mod_add(res, mod_mul(mod_norm(v), seg));
            l = r + 1;
        }
        return res;
    }

    int64 S(int64 n) {
        int64 ans = mod_add(mod_norm(n), T(n));
        return mod_mul(ans, INV2);
    }

    static int64 lcm_ll(int64 a, int64 b) {
        return a / std::gcd(a, b) * b;
    }

    static int64 brute_S(int n) {
        int64 total = 0;
        for (int k = 1; k <= n; ++k) {
            int64 sum = 0;
            for (int i = 1; i <= k; ++i) sum += lcm_ll(k, i);
            if (sum % k != 0) {
                std::cerr << "[VALIDATION] lcm-sum not divisible by k=" << k << "\n";
                std::exit(1);
            }
            total += sum / k;
        }
        return total;
    }

    void run_validations() {
        int64 brute_100 = brute_S(100);
        if (brute_100 != 122726LL) {
            std::cerr << "[VALIDATION] brute S(100)=" << brute_100 << " expected 122726\n";
            std::exit(1);
        }

        const int checks[] = {1, 2, 3, 10, 50, 100, 200};
        for (int n : checks) {
            int64 brute = brute_S(n) % MOD;
            int64 fast = S(n);
            if (brute != fast) {
                std::cerr << "[VALIDATION] mismatch S(" << n << "): brute=" << brute
                          << " fast=" << fast << "\n";
                std::exit(1);
            }
        }

        int max_check = std::min(LIM, 200000);
        int64 direct = 0;
        for (int i = 1; i <= max_check; ++i) {
            direct = mod_add(direct, mod_mul(static_cast<int64>(mu[i]), static_cast<int64>(i)));
            if (H(i) != direct) {
                std::cerr << "[VALIDATION] H(" << i << ") mismatch\n";
                std::exit(1);
            }
        }

        for (int i = 1; i <= max_check; i += 97) {
            if (G(i) != prefG[i]) {
                std::cerr << "[VALIDATION] G(" << i << ") mismatch\n";
                std::exit(1);
            }
        }
    }
};

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const int64 N = DEFAULT_N;
    Solver solver(N);
    solver.run_validations();
    std::cout << solver.S(N) << "\n";
    return 0;
}
