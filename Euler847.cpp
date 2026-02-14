#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <functional>
#include <vector>

using namespace std;

namespace {

static constexpr int64_t MOD = 1'000'000'007LL;

int64_t mod_add(int64_t a, int64_t b) {
    int64_t r = a + b;
    if (r >= MOD) r -= MOD;
    return r;
}

int64_t mod_sub(int64_t a, int64_t b) {
    int64_t r = a - b;
    if (r < 0) r += MOD;
    return r;
}

int64_t mod_mul(int64_t a, int64_t b) {
    return static_cast<int64_t>((__int128)a * b % MOD);
}

int64_t mod_pow(int64_t a, int64_t e) {
    int64_t r = 1 % MOD;
    int64_t x = (a % MOD + MOD) % MOD;
    while (e > 0) {
        if (e & 1) r = mod_mul(r, x);
        x = mod_mul(x, x);
        e >>= 1;
    }
    return r;
}

const int64_t INV2 = mod_pow(2, MOD - 2);
const int64_t INV6 = mod_pow(6, MOD - 2);

int64_t sum1_mod(long long n) {
    if (n < 0) return 0;
    int64_t n_mod = n % MOD;
    return mod_mul(mod_mul(n_mod, (n_mod + 1) % MOD), INV2);
}

int64_t sum2_mod(long long n) {
    if (n < 0) return 0;
    int64_t n_mod = n % MOD;
    int64_t n1 = (n_mod + 1) % MOD;
    int64_t two_n1 = (mod_mul(2 % MOD, n_mod) + 1) % MOD;
    return mod_mul(mod_mul(mod_mul(n_mod, n1), two_n1), INV6);
}

int64_t sum1_range_mod(long long l, long long r) {
    if (l > r) return 0;
    return mod_sub(sum1_mod(r), sum1_mod(l - 1));
}

int64_t sum2_range_mod(long long l, long long r) {
    if (l > r) return 0;
    return mod_sub(sum2_mod(r), sum2_mod(l - 1));
}

// Sum_{s=l..r} (s+1)(N-s+1).
int64_t sum_pair_mod(long long l, long long r, long long N) {
    if (l > r) return 0;
    long long cnt = r - l + 1;
    int64_t S1 = sum1_range_mod(l, r);
    int64_t S2 = sum2_range_mod(l, r);
    int64_t sum_s_plus = mod_add(S1, cnt % MOD);
    int64_t sum_s2_plus = mod_add(S2, S1);
    return mod_sub(mod_mul((N + 1) % MOD, sum_s_plus), sum_s2_plus);
}

// Sum_{s=l..r} (s+1)(s+2)/2.
int64_t sum_tri_mod(long long l, long long r) {
    if (l > r) return 0;
    long long cnt = r - l + 1;
    int64_t S1 = sum1_range_mod(l, r);
    int64_t S2 = sum2_range_mod(l, r);
    int64_t total = mod_add(S2, mod_add(mod_mul(3, S1), mod_mul(2, cnt % MOD)));
    return mod_mul(total, INV2);
}

// Sum_{a=l..r} (D-a+1)^2.
int64_t sum_M1_sq_mod(long long l, long long r, long long D) {
    if (l > r) return 0;
    long long u_low = D - r + 1;
    long long u_high = D - l + 1;
    return sum2_range_mod(u_low, u_high);
}

// Sum of t(t+1)/2 for arithmetic progression: t = t0 + step*i, i=0..n-1.
int64_t sum_t_t1_mod(long long t0, long long step, long long n) {
    if (n <= 0) return 0;
    int64_t n_mod = n % MOD;
    int64_t n1_mod = (n - 1) % MOD;
    if (n1_mod < 0) n1_mod += MOD;
    int64_t sum_i = mod_mul(mod_mul(n_mod, n1_mod), INV2);
    int64_t two_n1 = (mod_mul(2 % MOD, n1_mod) + 1) % MOD;
    int64_t sum_i2 = mod_mul(mod_mul(mod_mul(n_mod, n1_mod), two_n1), INV6);

    int64_t t0_mod = t0 % MOD;
    if (t0_mod < 0) t0_mod += MOD;
    int64_t step_mod = step % MOD;
    if (step_mod < 0) step_mod += MOD;

    int64_t sum_t = mod_add(mod_mul(n_mod, t0_mod), mod_mul(step_mod, sum_i));
    int64_t t0_sq = mod_mul(t0_mod, t0_mod);
    int64_t term1 = mod_mul(n_mod, t0_sq);
    int64_t term2 = mod_mul(mod_mul(mod_mul(2, t0_mod), step_mod), sum_i);
    int64_t term3 = mod_mul(mod_mul(step_mod, step_mod), sum_i2);
    int64_t sum_t2 = mod_add(term1, mod_add(term2, term3));

    return mod_mul(mod_add(sum_t2, sum_t), INV2);
}

int64_t T_mod(long long x) {
    if (x < 0) return 0;
    int64_t xm = x % MOD;
    return mod_mul(mod_mul(mod_mul((xm + 1) % MOD, (xm + 2) % MOD), (xm + 3) % MOD), INV6);
}

// Count triples 0<=a,b,c<=D, a+b+c<=N, and a+b>D, a+c>D, b+c>D.
int64_t count_within_cube_mod(long long N, long long D) {
    if (N < 0) return 0;
    if (N > 2 * D) N = 2 * D;

    int64_t U;
    if (N <= D) {
        U = T_mod(N);
    } else {
        U = mod_sub(T_mod(N), mod_mul(3, T_mod(N - D - 1)));
    }

    int64_t A;
    if (N <= D) {
        A = sum_pair_mod(0, N, N);
    } else {
        long long s1 = N - D;
        int64_t part1 = mod_mul((D + 1) % MOD,
                                mod_mul((s1 + 1) % MOD, (s1 + 2) % MOD) * INV2 % MOD);
        int64_t part2 = sum_pair_mod(s1 + 1, D, N);
        A = mod_add(part1, part2);
    }

    int64_t AB;
    if (N <= D) {
        AB = sum_tri_mod(0, N);
    } else {
        long long a0 = 2 * D - N;
        int64_t mid1 = sum_M1_sq_mod(0, a0 - 1, D);
        int64_t mid2 = sum_t_t1_mod(a0, -1, a0);
        int64_t mid = mod_sub(mid1, mid2);
        int64_t full = sum_M1_sq_mod(a0, D, D);
        AB = mod_add(mid, full);
    }

    int64_t ABC;
    if (N <= D) {
        ABC = sum_tri_mod(0, N);
    } else {
        long long A1 = N - D;
        long long a_half = D / 2;
        int64_t total = 0;

        long long r1_end = min(A1, a_half - 1);
        if (r1_end >= 0) {
            int64_t mid1 = sum_M1_sq_mod(0, r1_end, D);
            long long n = r1_end + 1;
            int64_t mid2 = sum_t_t1_mod(D, -2, n);
            total = mod_add(total, mod_sub(mid1, mid2));
        }

        long long r1_full_start = max(0LL, a_half);
        if (r1_full_start <= A1) {
            total = mod_add(total, sum_M1_sq_mod(r1_full_start, A1, D));
        }

        long long l2 = A1 + 1;
        if (l2 <= D) {
            long long a0 = 2 * D - N;
            long long r2_end = min(a0 - 1, D);
            if (l2 <= r2_end) {
                int64_t mid1 = sum_M1_sq_mod(l2, r2_end, D);
                long long n = r2_end - l2 + 1;
                int64_t mid2 = sum_t_t1_mod(a0 - l2, -1, n);
                total = mod_add(total, mod_sub(mid1, mid2));
            }

            long long l_full = max(l2, a0);
            if (l_full <= D) {
                total = mod_add(total, sum_M1_sq_mod(l_full, D, D));
            }
        }

        ABC = total;
    }

    int64_t res = mod_sub(mod_add(mod_sub(U, mod_mul(3, A)), mod_mul(3, AB)), ABC);
    return res;
}

vector<unordered_map<long long, int64_t>> memo;

int64_t hard_count_mod(int k, long long N) {
    if (N < 0) return 0;
    if (k == 1) {
        if (N <= 2) return 0;
        return mod_sub(T_mod(N), T_mod(2));
    }
    if (N > (1LL << k)) N = 1LL << k;
    auto& mp = memo[k];
    auto it = mp.find(N);
    if (it != mp.end()) return it->second;

    long long D = 1LL << (k - 1);
    int64_t res;
    if (N <= D) {
        res = count_within_cube_mod(N, D);
    } else {
        res = mod_add(count_within_cube_mod(N, D),
                      mod_mul(3, hard_count_mod(k - 1, N - D)));
    }
    mp.emplace(N, res);
    return res;
}

int64_t sum_C_range(long long L, long long R) {
    if (L > R) return 0;
    long long cnt = R - L + 1;
    int64_t S1 = sum1_range_mod(L, R);
    int64_t S2 = sum2_range_mod(L, R);
    int64_t total = mod_add(S2, mod_add(mod_mul(3, S1), mod_mul(2, cnt % MOD)));
    return mod_mul(total, INV2);
}

int64_t compute_H(long long N) {
    if (N <= 1) return 0;
    int max_k = 64 - __builtin_clzll(static_cast<unsigned long long>(N - 1));

    int64_t base = 0;
    for (int k = 1; k <= max_k; ++k) {
        long long L = (1LL << (k - 1)) + 1;
        long long R = min(N, 1LL << k);
        if (L > R) continue;
        int64_t sumC = sum_C_range(L, R);
        base = (base + mod_mul(k % MOD, sumC)) % MOD;
    }

    memo.assign(max_k + 1, {});
    int64_t extra = 0;
    for (int k = 1; k <= max_k; ++k) {
        long long X = min(N, 1LL << k);
        extra = mod_add(extra, hard_count_mod(k, X));
    }

    return mod_add(base, extra);
}

void run_validation() {
    struct Case { long long N; long long expected; };
    const Case cases[] = {
        {6, 203},
        {20, 7718},
        {111, 1634144},
    };
    for (const auto& c : cases) {
        int64_t got = compute_H(c.N);
        if (got != (c.expected % MOD)) {
            cerr << "Validation failed for N=" << c.N
                 << ": got " << got << ", expected " << c.expected << "\n";
            exit(1);
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long N = 0;
    for (int i = 0; i < 19; ++i) N = N * 10 + 1; // R_19
    bool validate = true;

    // Optional CLI: ./a.out [N] [validate(0/1)]
    if (argc >= 2) N = stoll(argv[1]);
    if (argc >= 3) validate = (stoi(argv[2]) != 0);

    if (validate) run_validation();

    cout << compute_H(N) << "\n";
    return 0;
}
