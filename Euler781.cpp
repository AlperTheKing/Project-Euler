#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using std::int64_t;

namespace {

constexpr int MOD = 1000000007;

inline int add_mod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

inline int sub_mod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

inline int mul_mod(int64_t a, int64_t b) {
    return static_cast<int>(a * b % MOD);
}

void log_series(const std::vector<int>& coeff,
                const std::vector<int>& inv,
                std::vector<int>& logc) {
    const int L = static_cast<int>(coeff.size()) - 1;
    logc.assign(L + 1, 0);
    for (int n = 1; n <= L; ++n) {
        int64_t s = 0;
        for (int k = 1; k < n; ++k) {
            s += static_cast<int64_t>(k) * logc[k] % MOD * coeff[n - k] % MOD;
            if (s >= MOD) s %= MOD;
        }
        s %= MOD;
        int64_t val = coeff[n] - s * inv[n] % MOD;
        if (val < 0) val += MOD;
        logc[n] = static_cast<int>(val);
    }
}

void build_coeff_even(int i, int maxm, const std::vector<int>& inv, std::vector<int>& coeff) {
    coeff.assign(maxm + 1, 0);
    coeff[0] = 1;
    int64_t S = 1;
    int64_t term = 1;
    int64_t M_prev2 = 1;  // M(0)
    int64_t M_prev1 = 1;  // M(1)
    const int inv_i = inv[i];
    for (int m = 1; m <= maxm; ++m) {
        term = term * (MOD - inv_i) % MOD;
        term = term * inv[m] % MOD;
        S = (S + term) % MOD;
        int64_t M = (M_prev1 + static_cast<int64_t>(m - 1) * i % MOD * M_prev2) % MOD;
        coeff[m] = mul_mod(S, M);
        M_prev2 = M_prev1;
        M_prev1 = M;
    }
}

void build_coeff_odd(int i, int maxr, const std::vector<int>& inv, std::vector<int>& coeff) {
    const int mlimit = 2 * maxr;
    coeff.assign(maxr + 1, 0);
    coeff[0] = 1;
    int64_t S = 1;
    int64_t term = 1;
    int64_t M_even_prev = 1;  // M(0)
    const int inv_i = inv[i];
    for (int m = 1; m <= mlimit; ++m) {
        term = term * (MOD - inv_i) % MOD;
        term = term * inv[m] % MOD;
        S = (S + term) % MOD;
        if ((m & 1) == 0) {
            int64_t M_even = static_cast<int64_t>(m - 1) * i % MOD * M_even_prev % MOD;
            const int r = m / 2;
            coeff[r] = mul_mod(S, M_even);
            M_even_prev = M_even;
        }
    }
}

void build_coeff_i1_star(int N,
                         const std::vector<int>& inv,
                         std::vector<int>& coeff,
                         std::vector<int>& coeff_star) {
    const int maxr = N;
    const int mlimit = 2 * maxr;
    coeff.assign(maxr + 1, 0);
    coeff_star.assign(maxr + 1, 0);
    coeff[0] = 1;
    coeff_star[0] = 1;
    int64_t S = 1;
    int64_t term = 1;
    int64_t S_prev = 0;
    int64_t M_even_prev = 1;
    for (int m = 1; m <= mlimit; ++m) {
        S_prev = S;
        term = term * (MOD - 1) % MOD;
        term = term * inv[m] % MOD;
        S = (S + term) % MOD;
        if ((m & 1) == 0) {
            int64_t M_even = static_cast<int64_t>(m - 1) * M_even_prev % MOD;
            const int r = m / 2;
            coeff[r] = mul_mod(S, M_even);
            int64_t B = (static_cast<int64_t>(m + 1) * S + S_prev) % MOD;
            coeff_star[r] = mul_mod(B, M_even);
            M_even_prev = M_even;
        }
    }
}

void accumulate_logs(int tid,
                     int threads,
                     int n,
                     int N,
                     const std::vector<int>& inv,
                     std::vector<int>& out) {
    std::vector<int> coeff;
    std::vector<int> logc;

    for (int i = 1 + tid; i <= n; i += threads) {
        if ((i & 1) == 0) {
            const int s = i / 2;
            const int maxm = N / s;
            if (maxm == 0) continue;
            build_coeff_even(i, maxm, inv, coeff);
            if (maxm == 1) {
                out[s] = add_mod(out[s], coeff[1]);
                continue;
            }
            log_series(coeff, inv, logc);
            for (int m = 1; m <= maxm; ++m) {
                const int idx = s * m;
                out[idx] = add_mod(out[idx], logc[m]);
            }
        } else {
            const int s = i;
            const int maxr = N / s;
            if (maxr == 0) continue;
            build_coeff_odd(i, maxr, inv, coeff);
            if (maxr == 1) {
                out[s] = add_mod(out[s], coeff[1]);
                continue;
            }
            log_series(coeff, inv, logc);
            for (int r = 1; r <= maxr; ++r) {
                const int idx = s * r;
                out[idx] = add_mod(out[idx], logc[r]);
            }
        }
    }
}

std::vector<int> exp_series(const std::vector<int>& logc, const std::vector<int>& inv) {
    const int N = static_cast<int>(logc.size()) - 1;
    std::vector<int> res(N + 1, 0);
    res[0] = 1;
    std::vector<int> kl(N + 1, 0);
    for (int k = 1; k <= N; ++k) {
        kl[k] = mul_mod(k, logc[k]);
    }
    for (int n = 1; n <= N; ++n) {
        int64_t s = 0;
        for (int k = 1; k <= n; ++k) {
            s += static_cast<int64_t>(kl[k]) * res[n - k] % MOD;
            if (s >= MOD) s %= MOD;
        }
        res[n] = mul_mod(s % MOD, inv[n]);
    }
    return res;
}

std::vector<int> divide_series(const std::vector<int>& A, const std::vector<int>& D) {
    const int N = static_cast<int>(A.size()) - 1;
    std::vector<int> C(N + 1, 0);
    C[0] = 1;
    for (int n = 1; n <= N; ++n) {
        int64_t s = A[n];
        for (int k = 1; k <= n; ++k) {
            s -= static_cast<int64_t>(C[n - k]) * D[k] % MOD;
            if (s < 0) s += MOD;
        }
        C[n] = static_cast<int>(s % MOD);
    }
    return C;
}

}  // namespace

int main() {
    const int n = 50000;
    if (n & 1) {
        std::cout << 0 << '\n';
        return 0;
    }
    const int N = n / 2;

    const int max_inv = 2 * N;
    std::vector<int> inv(max_inv + 1, 0);
    inv[1] = 1;
    for (int i = 2; i <= max_inv; ++i) {
        inv[i] = MOD - static_cast<int64_t>(MOD / i) * inv[MOD % i] % MOD;
    }

    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<int>> locals(threads, std::vector<int>(N + 1, 0));
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            accumulate_logs(static_cast<int>(t), static_cast<int>(threads), n, N, inv, locals[t]);
        });
    }
    for (auto &th : workers) th.join();

    std::vector<int> logD(N + 1, 0);
    for (int idx = 0; idx <= N; ++idx) {
        int64_t sum = 0;
        for (unsigned t = 0; t < threads; ++t) {
            sum += locals[t][idx];
        }
        logD[idx] = static_cast<int>(sum % MOD);
    }

    std::vector<int> D = exp_series(logD, inv);

    // Adjust log for i=1 to build logA.
    std::vector<int> coeff1;
    std::vector<int> coeff1_star;
    std::vector<int> log1;
    std::vector<int> log1_star;
    build_coeff_i1_star(N, inv, coeff1, coeff1_star);
    log_series(coeff1, inv, log1);
    log_series(coeff1_star, inv, log1_star);

    std::vector<int> logA = logD;
    for (int i = 1; i <= N; ++i) {
        logA[i] = add_mod(sub_mod(logA[i], log1[i]), log1_star[i]);
    }

    std::vector<int> A = exp_series(logA, inv);
    std::vector<int> C = divide_series(A, D);

    if (N >= 2 && C[2] != 5) {
        std::cerr << "Validation failed for F(4)." << '\n';
        return 1;
    }
    if (N >= 4 && C[4] != 319) {
        std::cerr << "Validation failed for F(8)." << '\n';
        return 1;
    }

    std::cout << C[N] << '\n';
    return 0;
}
