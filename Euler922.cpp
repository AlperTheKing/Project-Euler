#include <array>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using int64 = long long;

static constexpr int64 MOD = 1000000007LL;
static constexpr int XOR_SIZE = 64;

using Row = std::array<int64, XOR_SIZE>;

static inline int64 mul_mod(int64 a, int64 b) {
    return static_cast<int64>((static_cast<__int128>(a) * b) % MOD);
}

static void fwt_xor(Row& a, bool invert) {
    for (int len = 1; len < XOR_SIZE; len <<= 1) {
        for (int i = 0; i < XOR_SIZE; i += (len << 1)) {
            for (int j = 0; j < len; ++j) {
                int64 u = a[i + j];
                int64 v = a[i + j + len];
                int64 x = u + v;
                if (x >= MOD) {
                    x -= MOD;
                }
                int64 y = u - v;
                if (y < 0) {
                    y += MOD;
                }
                a[i + j] = x;
                a[i + j + len] = y;
            }
        }
    }
    if (!invert) {
        return;
    }
    const int64 inv2 = (MOD + 1) / 2;
    for (int len = 1; len < XOR_SIZE; len <<= 1) {
        for (int i = 0; i < XOR_SIZE; i += (len << 1)) {
            for (int j = 0; j < len; ++j) {
                a[i + j] = (a[i + j] * inv2) % MOD;
                a[i + j + len] = (a[i + j + len] * inv2) % MOD;
            }
        }
    }
}

static int64 compute_R(int m, int w) {
    const int maxD = w - 2;
    const int D = 2 * maxD + 1;
    const int S = 2 * maxD * m + 1;
    const int offset = maxD * m;

    std::vector<Row> f(D);
    for (int i = 0; i < D; ++i) {
        f[i].fill(0);
    }

    for (int a = 1; a <= w - 2; ++a) {
        for (int b = 1; b <= w - a - 1; ++b) {
            const int maxk = w - a - b;
            const int d_idx = b - a + maxD;
            Row& row = f[d_idx];
            for (int g = 0; g < maxk; ++g) {
                row[g] = (row[g] + 1) % MOD;
            }
        }
    }

    for (int d = 0; d < D; ++d) {
        fwt_xor(f[d], false);
    }

    std::vector<Row> dp(S), dp_hat(S), new_hat(S), new_dp(S);
    for (int i = 0; i < S; ++i) {
        dp[i].fill(0);
    }
    dp[offset][0] = 1;

    std::vector<int> delta(D);
    for (int d = 0; d < D; ++d) {
        delta[d] = d - maxD;
    }

    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads < 1) {
        threads = 1;
    }
    if (threads > XOR_SIZE) {
        threads = XOR_SIZE;
    }

    for (int step = 0; step < m; ++step) {
        dp_hat = dp;
        for (int s = 0; s < S; ++s) {
            fwt_xor(dp_hat[s], false);
        }

        for (int s = 0; s < S; ++s) {
            new_hat[s].fill(0);
        }

        auto worker = [&](int freq_start, int freq_end) {
            for (int freq = freq_start; freq < freq_end; ++freq) {
                for (int sd = 0; sd < S; ++sd) {
                    int64 acc = 0;
                    for (int d = 0; d < D; ++d) {
                        const int sd_prev = sd - delta[d];
                        if (sd_prev < 0 || sd_prev >= S) {
                            continue;
                        }
                        acc += mul_mod(dp_hat[sd_prev][freq], f[d][freq]);
                        if (acc >= MOD) {
                            acc -= MOD;
                        }
                    }
                    new_hat[sd][freq] = acc;
                }
            }
        };

        if (threads == 1) {
            worker(0, XOR_SIZE);
        } else {
            std::vector<std::thread> pool;
            pool.reserve(threads);
            const int chunk = (XOR_SIZE + threads - 1) / threads;
            for (int t = 0; t < threads; ++t) {
                const int start = t * chunk;
                const int end = std::min(XOR_SIZE, start + chunk);
                if (start >= end) {
                    break;
                }
                pool.emplace_back(worker, start, end);
            }
            for (auto& th : pool) {
                th.join();
            }
        }

        new_dp = new_hat;
        for (int s = 0; s < S; ++s) {
            fwt_xor(new_dp[s], true);
        }
        dp.swap(new_dp);
    }

    int64 ans = 0;
    for (int sd = 0; sd < S; ++sd) {
        const int sum_d = sd - offset;
        int64 row_sum = 0;
        for (int x = 0; x < XOR_SIZE; ++x) {
            row_sum += dp[sd][x];
            if (row_sum >= MOD) {
                row_sum -= MOD;
            }
        }
        if (sum_d > 0) {
            ans += row_sum;
        } else if (sum_d == 0) {
            int64 zero = dp[sd][0];
            int64 add = row_sum - zero;
            if (add < 0) {
                add += MOD;
            }
            ans += add;
        }
        if (ans >= MOD) {
            ans %= MOD;
        }
    }
    return ans % MOD;
}

static bool run_validation() {
    struct Case {
        int m;
        int w;
        int64 expected;
    };
    const Case cases[] = {
        {2, 4, 7},
        {3, 9, 314104},
    };
    for (const auto& tc : cases) {
        int64 got = compute_R(tc.m, tc.w);
        if (got != tc.expected) {
            std::cerr << "Validation failed for R(" << tc.m << "," << tc.w << ") expected="
                      << tc.expected << " got=" << got << "\n";
            return false;
        }
    }
    return true;
}

int main() {
    if (!run_validation()) {
        return 1;
    }
    const int m = 8;
    const int w = 64;
    std::cout << compute_R(m, w) << "\n";
    return 0;
}
