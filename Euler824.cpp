#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using u64 = unsigned long long;
using u128 = unsigned __int128;

namespace {

constexpr u64 P = 10000019ULL;
constexpr u64 MOD = P * P;

std::vector<u64> fact;
std::vector<u64> H;
std::vector<u64> inv_small;
u64 fact_p_minus_1 = 0;

inline u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>(static_cast<u128>(a) * b % MOD);
}

u64 pow_mod(u64 base, u64 exp) {
    u64 res = 1;
    while (exp > 0) {
        if (exp & 1ULL) res = mul_mod(res, base);
        base = mul_mod(base, base);
        exp >>= 1ULL;
    }
    return res;
}

u64 modinv(u64 a, u64 mod) {
    long long t = 0, new_t = 1;
    long long r = static_cast<long long>(mod);
    long long new_r = static_cast<long long>(a);
    while (new_r != 0) {
        long long q = r / new_r;
        long long tmp_t = static_cast<long long>(static_cast<__int128>(t) - static_cast<__int128>(q) * new_t);
        long long tmp_r = r - q * new_r;
        t = new_t;
        new_t = tmp_t;
        r = new_r;
        new_r = tmp_r;
    }
    if (r != 1) return 0;
    if (t < 0) t += static_cast<long long>(mod);
    return static_cast<u64>(t);
}

u64 vp_factorial(u64 n) {
    u64 e = 0;
    while (n > 0) {
        n /= P;
        e += n;
    }
    return e;
}

u64 fact_without_p(u64 n) {
    if (n == 0) return 1;
    u64 q = n / P;
    u64 r = n % P;
    u64 res = fact_without_p(q);
    res = mul_mod(res, pow_mod(fact_p_minus_1, q));
    res = mul_mod(res, fact[static_cast<size_t>(r)]);
    if (r != 0) {
        u64 q_mod = q % P;
        u64 factor = 1;
        if (q_mod != 0) {
            u64 tmp = mul_mod(mul_mod(q_mod, P), H[static_cast<size_t>(r)]);
            factor = (1 + tmp) % MOD;
        }
        res = mul_mod(res, factor);
    }
    return res;
}

u64 binom_mod_p2(u64 n, u64 k) {
    if (k > n) return 0;
    u64 e = vp_factorial(n);
    e -= vp_factorial(k);
    e -= vp_factorial(n - k);
    if (e >= 2) return 0;
    u64 a = fact_without_p(n);
    u64 b = fact_without_p(k);
    u64 c = fact_without_p(n - k);
    u64 denom = mul_mod(b, c);
    u64 inv_denom = modinv(denom, MOD);
    u64 res = mul_mod(a, inv_denom);
    if (e == 1) res = mul_mod(res, P);
    return res;
}

u64 coeff_lambda_plus(u64 m, u64 k) {
    if (k == 0) return 1;
    if (m == 0) return 0;
    // Lagrange inversion: [x^k] lambda_+^m = (m/k) * C(m-k-1, k-1).
    u64 n = m - k - 1;
    u64 r = k - 1;
    u64 c = binom_mod_p2(n, r);
    u64 inv_k = modinv(k % MOD, MOD);
    return mul_mod(c, mul_mod(m % MOD, inv_k));
}

void init_tables(u64 max_j) {
    fact.assign(P, 0);
    fact[0] = 1;
    for (u64 i = 1; i < P; ++i) {
        fact[static_cast<size_t>(i)] = mul_mod(fact[static_cast<size_t>(i - 1)], i);
    }
    fact_p_minus_1 = fact[static_cast<size_t>(P - 1)];

    std::vector<u64> inv(P, 0);
    inv[1] = 1;
    for (u64 i = 2; i < P; ++i) {
        inv[static_cast<size_t>(i)] =
            (MOD - mul_mod(MOD / i, inv[static_cast<size_t>(MOD % i)])) % MOD;
    }

    H.assign(P, 0);
    for (u64 i = 1; i < P; ++i) {
        H[static_cast<size_t>(i)] = (H[static_cast<size_t>(i - 1)] + inv[static_cast<size_t>(i)]) % MOD;
    }

    inv_small.assign(static_cast<size_t>(max_j + 1), 0);
    if (max_j >= 1) inv_small[1] = inv[1];
    for (u64 i = 2; i <= max_j; ++i) {
        inv_small[static_cast<size_t>(i)] = inv[static_cast<size_t>(i)];
    }
}

u64 compute_L(u64 N, u64 K) {
    u64 J = K / N;
    if (J >= inv_small.size()) {
        std::cerr << "Precomputed inverses are too small." << '\n';
        std::exit(1);
    }

    std::vector<u64> binomN(static_cast<size_t>(J + 1), 0);
    binomN[0] = 1;
    for (u64 j = 0; j < J; ++j) {
        u64 num = (N - j) % MOD;
        u64 next = mul_mod(binomN[static_cast<size_t>(j)], mul_mod(num, inv_small[static_cast<size_t>(j + 1)]));
        binomN[static_cast<size_t>(j + 1)] = next;
    }

    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    if (threads > J + 1) threads = static_cast<unsigned>(J + 1);
    u64 chunk = (J + 1 + threads - 1) / threads;

    std::vector<u64> partial(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        u64 start = static_cast<u64>(t) * chunk;
        u64 end = std::min<u64>(J + 1, start + chunk);
        workers.emplace_back([&, t, start, end]() {
            u64 sum = 0;
            for (u64 j = start; j < end; ++j) {
                u64 k = K - N * j;
                u64 coeff = coeff_lambda_plus(N * (N - 2 * j), k);
                u64 term = mul_mod(binomN[static_cast<size_t>(j)], coeff);
                if ((N & 1ULL) && (j & 1ULL) && term != 0) {
                    term = MOD - term;
                }
                sum += term;
                if (sum >= MOD) sum -= MOD;
            }
            partial[t] = sum;
        });
    }

    for (auto &th : workers) th.join();

    u64 total = 0;
    for (u64 val : partial) {
        total += val;
        if (total >= MOD) total -= MOD;
    }
    return total;
}

void validate() {
    const u64 check1 = compute_L(2, 2);
    if (check1 != 4) {
        std::cerr << "Validation failed for L(2,2)." << '\n';
        std::exit(1);
    }
    const u64 check2 = compute_L(6, 12);
    if (check2 != 4204761ULL) {
        std::cerr << "Validation failed for L(6,12)." << '\n';
        std::exit(1);
    }
}

}  // namespace

int main() {
    const u64 N = 1000000000ULL;
    const u64 K = 1000000000000000ULL;
    const u64 max_j = K / N;

    init_tables(max_j);
    validate();

    u64 result = compute_L(N, K);
    std::cout << result << '\n';
    return 0;
}
