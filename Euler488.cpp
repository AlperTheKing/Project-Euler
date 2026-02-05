#include <array>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr uint64_t MOD = 1'000'000'000ULL;

inline uint64_t add_mod(uint64_t a, uint64_t b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}

inline uint64_t sub_mod(uint64_t a, uint64_t b) {
    return (a >= b) ? (a - b) : (a + MOD - b);
}

inline uint64_t mul_mod(uint64_t a, uint64_t b) {
    return static_cast<uint64_t>((__int128)(a % MOD) * (b % MOD) % MOD);
}

uint64_t sum_linear_mod(uint64_t m) {
    // 1 + 2 + ... + m
    if (m == 0) return 0;
    uint64_t a = m, b = m + 1;
    if ((a & 1ULL) == 0) {
        a >>= 1;
    } else {
        b >>= 1;
    }
    return mul_mod(a, b);
}

uint64_t sum_square_mod(uint64_t m) {
    // 1^2 + 2^2 + ... + m^2 = m(m+1)(2m+1)/6
    if (m == 0) return 0;
    uint64_t a = m, b = m + 1, c = 2 * m + 1;

    // divide by 2
    if ((a & 1ULL) == 0) {
        a >>= 1;
    } else if ((b & 1ULL) == 0) {
        b >>= 1;
    } else {
        c >>= 1;
    }

    // divide by 3
    if (a % 3ULL == 0) {
        a /= 3ULL;
    } else if (b % 3ULL == 0) {
        b /= 3ULL;
    } else {
        c /= 3ULL;
    }

    return mul_mod(mul_mod(a, b), c);
}

uint64_t sum_choose2_mod(uint64_t m) {
    // sum_{n=1..m} n(n-1)/2 = m(m+1)(m-1)/6
    if (m <= 1) return 0;
    uint64_t a = m, b = m + 1, c = m - 1;

    // divide by 2
    if ((a & 1ULL) == 0) {
        a >>= 1;
    } else if ((b & 1ULL) == 0) {
        b >>= 1;
    } else {
        c >>= 1;
    }

    // divide by 3
    if (a % 3ULL == 0) {
        a /= 3ULL;
    } else if (b % 3ULL == 0) {
        b /= 3ULL;
    } else {
        c /= 3ULL;
    }

    return mul_mod(mul_mod(a, b), c);
}

struct Block {
    uint64_t p;
    uint64_t m;
};

class Solver488 {
public:
    Solver488() {
        precompute_prefix_xor_sums();
    }

    uint64_t solve_mod(uint64_t N, bool allow_multithreading = true) const {
        vector<Block> blocks;
        blocks.reserve(64);
        for (int k = 1; k <= 62; ++k) {
            uint64_t p = (1ULL << k) - 1ULL;
            if (p + 1ULL >= N) break;
            uint64_t m = min<uint64_t>(p, N - 1ULL - p);
            blocks.push_back({p, m});
        }

        if (blocks.empty()) return 0;

        unsigned workers = 1;
        if (allow_multithreading) {
            unsigned hw = thread::hardware_concurrency();
            if (hw == 0) hw = 1;
            workers = min<unsigned>(hw, static_cast<unsigned>(blocks.size()));
            if (blocks.size() < 8) workers = 1;
        }

        if (workers == 1) {
            uint64_t ans = 0;
            for (const auto &blk : blocks) {
                ans = add_mod(ans, block_sum_mod(blk.p, blk.m));
            }
            return ans;
        }

        vector<uint64_t> partial(workers, 0);
        vector<thread> ts;
        ts.reserve(workers);
        for (unsigned tid = 0; tid < workers; ++tid) {
            ts.emplace_back([&, tid]() {
                size_t L = blocks.size() * tid / workers;
                size_t R = blocks.size() * (tid + 1) / workers;
                uint64_t local = 0;
                for (size_t i = L; i < R; ++i) {
                    local = add_mod(local, block_sum_mod(blocks[i].p, blocks[i].m));
                }
                partial[tid] = local;
            });
        }
        for (auto &th : ts) th.join();

        uint64_t ans = 0;
        for (uint64_t x : partial) ans = add_mod(ans, x);
        return ans;
    }

private:
    // S(n) = sum_{i=0..n-1} X(i),  X(i)=sum_{t=0..i-1}(i xor t)
    array<uint64_t, 63> s_pow2_mod{};
    array<uint64_t, 63> c_mod{};

    void precompute_prefix_xor_sums() {
        s_pow2_mod.fill(0);
        c_mod.fill(0);
        for (int m = 0; m <= 62; ++m) {
            uint64_t p = (1ULL << m);
            uint64_t x = p;
            uint64_t y = 3ULL * p - 1ULL;
            if ((x & 1ULL) == 0) {
                x >>= 1;
            } else {
                y >>= 1;
            }
            c_mod[m] = mul_mod(x, y); // C_m = p(3p-1)/2
        }
        for (int m = 0; m < 62; ++m) {
            uint64_t p = (1ULL << m) % MOD;
            uint64_t term = mul_mod(p, c_mod[m]);
            s_pow2_mod[m + 1] = add_mod(add_mod(s_pow2_mod[m], s_pow2_mod[m]), term);
        }
    }

    uint64_t prefix_sum_x_mod(uint64_t n) const {
        // returns S(n)=sum_{i=0..n-1}X(i)
        if (n <= 1) return 0;
        int msb = 63 - __builtin_clzll(n);
        uint64_t p = (1ULL << msb);
        if (n == p) return s_pow2_mod[msb];
        uint64_t r = n - p;
        uint64_t term = mul_mod(r % MOD, c_mod[msb]);
        return add_mod(add_mod(s_pow2_mod[msb], term), prefix_sum_x_mod(r));
    }

    uint64_t block_sum_mod(uint64_t p, uint64_t m) const {
        // P-positions in block p=2^k-1 are parameterized by:
        // c=p+n, b=p+t, a=(n xor t)-1 with 1<=n<=m and 0<=t<n.
        // Sum is then reduced to closed forms below.
        uint64_t s_x = prefix_sum_x_mod(m + 1); // sum_{n=1..m} X(n)
        uint64_t s1 = sum_linear_mod(m);        // sum n
        uint64_t s2 = sum_square_mod(m);        // sum n^2
        uint64_t sC2 = sum_choose2_mod(m);      // sum n(n-1)/2

        uint64_t two_p_minus_one = (2ULL * p - 1ULL) % MOD;

        uint64_t total = 0;
        total = add_mod(total, s_x);
        total = add_mod(total, s2);
        total = add_mod(total, sC2);
        total = add_mod(total, mul_mod(two_p_minus_one, s1));

        // Remove n odd term where t=n xor 1 = n-1 -> a=0 (not counted in F(N)).
        uint64_t odd_cnt = (m + 1ULL) >> 1;
        uint64_t odd_sq = mul_mod(odd_cnt % MOD, odd_cnt % MOD);
        uint64_t correction = add_mod(mul_mod(odd_cnt % MOD, two_p_minus_one),
                                      mul_mod(2ULL, odd_sq));
        return sub_mod(total, correction);
    }
};

uint64_t brute_F(uint32_t N) {
    auto idx = [N](uint32_t a, uint32_t b, uint32_t c) {
        return (a * N + b) * N + c;
    };

    vector<uint8_t> win(static_cast<size_t>(N) * N * N, 0);
    uint64_t ans = 0;

    for (uint32_t c = 2; c < N; ++c) {
        for (uint32_t b = 1; b < c; ++b) {
            for (uint32_t a = 0; a < b; ++a) {
                bool w = false;

                // Decrease a.
                for (uint32_t x = 0; x < a && !w; ++x) {
                    if (!win[idx(x, b, c)]) w = true;
                }

                // Decrease b.
                for (uint32_t x = 0; x < b && !w; ++x) {
                    if (x == a || x == c) continue;
                    uint32_t u, v, z;
                    if (x < a) {
                        u = x; v = a; z = c;
                    } else {
                        u = a; v = x; z = c;
                    }
                    if (!win[idx(u, v, z)]) w = true;
                }

                // Decrease c.
                for (uint32_t x = 0; x < c && !w; ++x) {
                    if (x == a || x == b) continue;
                    uint32_t u, v, z;
                    if (x < b) {
                        if (x < a) {
                            u = x; v = a; z = b;
                        } else {
                            u = a; v = x; z = b;
                        }
                    } else {
                        u = a; v = b; z = x;
                    }
                    if (!win[idx(u, v, z)]) w = true;
                }

                win[idx(a, b, c)] = static_cast<uint8_t>(w);
                if (!w && a > 0) {
                    ans += static_cast<uint64_t>(a) + b + c;
                }
            }
        }
    }
    return ans;
}

void validate() {
    Solver488 solver;

    const uint64_t f8 = solver.solve_mod(8, false);
    if (f8 != 42ULL) {
        cerr << "Validation failed: F(8)=" << f8 << " (expected 42)\n";
        exit(1);
    }

    const uint64_t f128 = solver.solve_mod(128, false);
    if (f128 != 496062ULL) {
        cerr << "Validation failed: F(128)=" << f128 << " (expected 496062)\n";
        exit(1);
    }

    // Extra structural checks against brute force.
    for (uint32_t n : {16u, 24u, 32u, 40u}) {
        uint64_t exact_brute = brute_F(n);
        uint64_t fast = solver.solve_mod(n, false);
        if (exact_brute != fast) {
            cerr << "Brute-force mismatch at N=" << n
                 << ": brute=" << exact_brute
                 << ", fast=" << fast << "\n";
            exit(1);
        }
    }
}

} // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    validate();

    Solver488 solver;
    const uint64_t ans = solver.solve_mod(1'000'000'000'000'000'000ULL, true);
    cout << setw(9) << setfill('0') << ans << "\n";
    return 0;
}
