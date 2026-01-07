#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

using std::int64_t;

namespace {

using i128 = __int128_t;

constexpr int64_t kTarget = 123567101113LL;
// Balance point between DFS (d <= B) and Pell (d > B). Must be > sqrt(n).
constexpr int64_t kBGlobal = 60000000LL;
constexpr int kMaxRoots = 128; // Enough for B=6e7 (max omega is 6 -> 64 roots).

struct PrimeRoot {
    int p;
    int64_t p2;
    int64_t r1;
    int64_t r2;
};

struct Roots {
    int size = 0;
    int64_t vals[kMaxRoots]{};
};

std::vector<PrimeRoot> g_primes;
std::vector<int> g_small_primes;

int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t res = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = (int64_t)((i128)res * base % mod);
        base = (int64_t)((i128)base * base % mod);
        exp >>= 1;
    }
    return res;
}

int64_t egcd(int64_t a, int64_t b, int64_t &x, int64_t &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int64_t x1 = 0, y1 = 0;
    int64_t g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (int64_t)((i128)y1 * (a / b));
    return g;
}

int64_t mod_inv(int64_t a, int64_t mod) {
    int64_t x = 0, y = 0;
    int64_t g = egcd(a, mod, x, y);
    if (g != 1) return 0;
    x %= mod;
    if (x < 0) x += mod;
    return x;
}

int64_t sqrt_minus_one_mod_p(int64_t p) {
    // Find any quadratic non-residue b, then b^((p-1)/4) is a sqrt(-1).
    for (int64_t b = 2; b < p; ++b) {
        if (mod_pow(b, (p - 1) / 2, p) == p - 1) {
            return mod_pow(b, (p - 1) / 4, p);
        }
    }
    return -1;
}

PrimeRoot build_prime_root(int p) {
    int64_t r0 = sqrt_minus_one_mod_p(p);
    int64_t p2 = (int64_t)p * p;
    int64_t m = (r0 * r0 + 1) / p;
    int64_t inv = mod_pow((2 * r0) % p, p - 2, p);
    int64_t k = (p - (m % p)) % p;
    k = (int64_t)((i128)k * inv % p);
    int64_t r = r0 + (int64_t)p * k;
    return PrimeRoot{p, p2, r, p2 - r};
}

void init_prime_roots(int64_t limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int64_t i = 2; i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (int64_t j = i * i; j <= limit; j += i) {
                is_prime[(size_t)j] = false;
            }
        }
    }
    for (int p = 5; p <= limit; ++p) {
        if (is_prime[p] && (p % 4 == 1)) {
            g_primes.push_back(build_prime_root(p));
        }
    }
}

void init_small_primes(int limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j <= limit; j += i) {
                is_prime[j] = false;
            }
        }
    }
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[i]) g_small_primes.push_back(i);
    }
}

int64_t count_solutions(int64_t n, int64_t mod, const Roots &roots) {
    int64_t q = n / mod;
    int64_t s = n - q * mod;
    int64_t count = q * roots.size;
    for (int i = 0; i < roots.size; ++i) {
        if (roots.vals[i] <= s) ++count;
    }
    return count;
}

void combine_roots(const Roots &roots, const PrimeRoot &pr, int64_t mod_old, Roots &out) {
    int64_t mod_new = pr.p2;
    int64_t inv = mod_inv(mod_old % mod_new, mod_new);
    out.size = roots.size * 2;
    if (out.size > kMaxRoots) {
        std::cerr << "Root buffer too small.\n";
        std::exit(1);
    }
    int idx = 0;
    for (int i = 0; i < roots.size; ++i) {
        int64_t ro = roots.vals[i];
        int64_t ro_mod = ro % mod_new;
        int64_t r_list[2] = {pr.r1, pr.r2};
        for (int j = 0; j < 2; ++j) {
            int64_t diff = r_list[j] - ro_mod;
            if (diff < 0) diff += mod_new;
            int64_t t = (int64_t)((i128)diff * inv % mod_new);
            int64_t new_r = ro + (int64_t)((i128)mod_old * t);
            out.vals[idx++] = new_r;
        }
    }
}

void dfs_part1(int start_idx,
               int64_t d,
               int64_t d2,
               const Roots &roots,
               int sign,
               int64_t n,
               int64_t B,
               int64_t &sum) {
    for (int i = start_idx; i < (int)g_primes.size(); ++i) {
        const PrimeRoot &pr = g_primes[i];
        if (pr.p > B) break;
        if (d > B / pr.p) break;
        int64_t next_d = d * pr.p;
        int64_t next_d2 = next_d * next_d;
        Roots next_roots;
        combine_roots(roots, pr, d2, next_roots);
        int64_t term = count_solutions(n, next_d2, next_roots);
        int next_sign = -sign;
        sum += (int64_t)next_sign * term;
        dfs_part1(i + 1, next_d, next_d2, next_roots, next_sign, n, B, sum);
    }
}

int64_t count_part1(int64_t n, int64_t B) {
    if (B < 5) return 0;
    int threads = (int)std::thread::hardware_concurrency();
    if (threads <= 0) threads = 4;
    std::vector<int64_t> sums(threads, 0);
    std::vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            int64_t local = 0;
            for (int idx = t; idx < (int)g_primes.size(); idx += threads) {
                const PrimeRoot &pr = g_primes[idx];
                if (pr.p > B) break;
                int64_t d = pr.p;
                int64_t d2 = d * d;
                Roots roots;
                roots.size = 2;
                roots.vals[0] = pr.r1;
                roots.vals[1] = pr.r2;
                int64_t term = count_solutions(n, d2, roots);
                local -= term; // mu(p) = -1
                dfs_part1(idx + 1, d, d2, roots, -1, n, B, local);
            }
            sums[t] = local;
        });
    }
    for (auto &th : workers) th.join();
    int64_t total = 0;
    for (int64_t v : sums) total += v;
    return total;
}

bool solve_negative_pell(int64_t k, int64_t n, int64_t &x, int64_t &y) {
    int64_t a0 = (int64_t)std::sqrt((long double)k);
    while ((a0 + 1) * (a0 + 1) <= k) ++a0;
    while (a0 * a0 > k) --a0;
    if (a0 * a0 == k) return false;

    int64_t m = 0, d = 1, a = a0;
    i128 p_prev1 = 1;
    i128 q_prev1 = 0;
    i128 p_curr = a0, q_curr = 1;
    if (p_curr > n || q_curr > n) return false;

    for (int step = 1; ; ++step) {
        m = d * a - m;
        d = (k - m * m) / d;
        a = (a0 + m) / d;

        i128 p_next = (i128)a * p_curr + p_prev1;
        i128 q_next = (i128)a * q_curr + q_prev1;

        p_prev1 = p_curr;
        q_prev1 = q_curr;
        p_curr = p_next;
        q_curr = q_next;

        if (d == 1 && a == 2 * a0) {
            if (step % 2 == 1 && p_prev1 <= n && q_prev1 <= n) {
                x = (int64_t)p_prev1;
                y = (int64_t)q_prev1;
                return true;
            }
            return false;
        }
        if (p_curr > n || q_curr > n) return false;
    }
}

int mobius_squarefree(int64_t y) {
    int mu = 1;
    int64_t tmp = y;
    for (int p : g_small_primes) {
        int64_t pp = (int64_t)p * p;
        if (pp > tmp) break;
        if (tmp % p == 0) {
            tmp /= p;
            if (tmp % p == 0) return 0;
            mu = -mu;
        }
    }
    if (tmp > 1) mu = -mu;
    return mu;
}

std::vector<char> build_valid_k(int64_t k_max) {
    std::vector<char> valid(k_max + 1, 1);
    if (k_max >= 0) valid[0] = 0;
    std::vector<bool> is_prime(k_max + 1, true);
    if (k_max >= 0) is_prime[0] = false;
    if (k_max >= 1) is_prime[1] = false;
    for (int64_t i = 2; i * i <= k_max; ++i) {
        if (is_prime[i]) {
            for (int64_t j = i * i; j <= k_max; j += i) is_prime[(size_t)j] = false;
        }
    }
    for (int64_t p = 2; p <= k_max; ++p) {
        if (!is_prime[p]) continue;
        if (p % 4 == 3) {
            for (int64_t j = p; j <= k_max; j += p) valid[(size_t)j] = 0;
        }
    }
    return valid;
}

int64_t count_part2(int64_t n, int64_t B) {
    int64_t k_max = (int64_t)(((i128)n * n + 1) / ((i128)B * B));
    if (k_max <= 0) return 0;
    std::vector<char> valid = build_valid_k(k_max);
    int threads = (int)std::thread::hardware_concurrency();
    if (threads <= 0) threads = 4;
    std::vector<int64_t> sums(threads, 0);
    std::vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            int64_t local = 0;
            for (int64_t k = 1 + t; k <= k_max; k += threads) {
                if (!valid[(size_t)k]) continue;
                int64_t s = (int64_t)std::sqrt((long double)k);
                if (s * s == k) continue;

                int64_t x0 = 0, y0 = 0;
                if (!solve_negative_pell(k, n, x0, y0)) continue;

                i128 X_mul = (i128)x0 * x0 + (i128)k * y0 * y0;
                i128 Y_mul = 2 * (i128)x0 * y0;
                i128 Y_mul_k = Y_mul * k;

                int64_t x = x0;
                int64_t y = y0;
                while (x <= n) {
                    if (y > B) {
                        int mu = mobius_squarefree(y);
                        if (mu != 0) local += mu;
                    }
                    if (X_mul > n || Y_mul_k > n) break;
                    i128 x_next = (i128)x * X_mul + (i128)y * Y_mul_k;
                    i128 y_next = (i128)x * Y_mul + (i128)y * X_mul;
                    if (x_next > n) break;
                    x = (int64_t)x_next;
                    y = (int64_t)y_next;
                }
            }
            sums[t] = local;
        });
    }
    for (auto &th : workers) th.join();
    int64_t total = 0;
    for (int64_t v : sums) total += v;
    return total;
}

int64_t compute_C(int64_t n) {
    int64_t B = std::min(kBGlobal, n);
    int64_t total = n; // d=1
    total += count_part1(n, B);
    total += count_part2(n, B);
    return total;
}

int64_t brute_force_C(int64_t n) {
    int64_t cnt = 0;
    for (int64_t x = 1; x <= n; ++x) {
        int64_t v = x * x + 1;
        bool squarefree = true;
        for (int p : g_small_primes) {
            int64_t pp = (int64_t)p * p;
            if (pp > v) break;
            if (v % pp == 0) {
                squarefree = false;
                break;
            }
        }
        if (squarefree) ++cnt;
    }
    return cnt;
}

} // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::cout << "Building prime roots...\n";
    init_prime_roots(kBGlobal);
    init_small_primes(400000);

    std::cout << "Validation checks...\n";
    int64_t c10 = compute_C(10);
    int64_t c10_bf = brute_force_C(10);
    std::cout << "C(10) = " << c10 << " (expected 9, brute " << c10_bf << ")\n";
    if (c10 != 9 || c10_bf != 9) {
        std::cerr << "Validation failed for C(10).\n";
        return 1;
    }

    int64_t c1000 = compute_C(1000);
    int64_t c1000_bf = brute_force_C(1000);
    std::cout << "C(1000) = " << c1000 << " (expected 895, brute " << c1000_bf << ")\n";
    if (c1000 != 895 || c1000_bf != 895) {
        std::cerr << "Validation failed for C(1000).\n";
        return 1;
    }

    std::cout << "Computing target...\n";
    int64_t result = compute_C(kTarget);
    std::cout << "C(" << kTarget << ") = " << result << "\n";
    return 0;
}
