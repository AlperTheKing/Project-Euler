#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr int kMod = 100000007;
constexpr int kFibLimit = 90;
constexpr int kDefaultBlockSize = 8192;

int mod_pow(int base, u64 exp, int mod) {
    u64 result = 1;
    u64 cur = static_cast<u64>(base % mod);
    while (exp > 0) {
        if (exp & 1ULL) result = (result * cur) % mod;
        cur = (cur * cur) % mod;
        exp >>= 1ULL;
    }
    return static_cast<int>(result);
}

struct BlockFactorial {
    int mod;
    int block_size;
    std::vector<int> block_fact;  // block_fact[b] = (b*B)! mod p (clipped to p-1).

    BlockFactorial(int mod_, int block_size_) : mod(mod_), block_size(block_size_) {
        if (block_size <= 0) block_size = kDefaultBlockSize;
        const int max_n = mod - 1;
        const int num_blocks = (max_n + block_size - 1) / block_size;
        block_fact.assign(num_blocks + 1, 1);

        for (int b = 1; b <= num_blocks; ++b) {
            const int start = (b - 1) * block_size + 1;
            const int end = std::min(b * block_size, max_n);
            u64 cur = static_cast<u64>(block_fact[b - 1]);
            for (int x = start; x <= end; ++x) {
                cur = (cur * static_cast<u64>(x)) % mod;
            }
            block_fact[b] = static_cast<int>(cur);
        }
    }

    int factorial_mod(int n) const {
        if (n < 0 || n >= mod) return 0;
        if (n == 0) return 1;

        const int b = n / block_size;
        u64 cur = static_cast<u64>(block_fact[b]);
        const int start = b * block_size + 1;
        for (int x = start; x <= n; ++x) {
            cur = (cur * static_cast<u64>(x)) % mod;
        }
        return static_cast<int>(cur);
    }

    int binom_small(int n, int k) const {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;

        const int nf = factorial_mod(n);
        const int kf = factorial_mod(k);
        const int nkf = factorial_mod(n - k);
        const int den = static_cast<int>((static_cast<u64>(kf) * nkf) % mod);
        const int inv_den = mod_pow(den, static_cast<u64>(mod - 2), mod);
        return static_cast<int>((static_cast<u64>(nf) * inv_den) % mod);
    }

    int binom_lucas(u64 n, u64 k) const {
        if (k > n) return 0;
        u64 nn = n;
        u64 kk = k;
        int result = 1;
        while (nn > 0 || kk > 0) {
            const int ni = static_cast<int>(nn % static_cast<u64>(mod));
            const int ki = static_cast<int>(kk % static_cast<u64>(mod));
            if (ki > ni) return 0;
            result = static_cast<int>((static_cast<u64>(result) * binom_small(ni, ki)) % mod);
            nn /= static_cast<u64>(mod);
            kk /= static_cast<u64>(mod);
        }
        return result;
    }

    int central_binom(u64 n) const { return binom_lucas(2ULL * n, n); }
};

bool attacks(int left_state, int right_state, int di, int dj) {
    static constexpr std::array<int, 4> a = {0, 0, 1, 1};
    static constexpr std::array<int, 4> b = {0, 1, 0, 1};

    const int dx = 2 * di + (a[right_state] - a[left_state]);
    const int dy = 2 * dj + (b[right_state] - b[left_state]);
    const int adx = std::abs(dx);
    const int ady = std::abs(dy);

    if (std::max(adx, ady) == 1) return true;                  // king move
    if ((adx == 1 && ady == 2) || (adx == 2 && ady == 1)) return true;  // knight move
    return false;
}

u64 exact_count_transfer(int n) {
    std::vector<std::vector<std::uint8_t>> rows;
    rows.reserve(static_cast<std::size_t>(4 * n));

    // Horizontal constraints imply this 4n-parameterized row family.
    for (int k = 0; k <= n; ++k) {
        if (k == 0) {
            for (int q = 0; q <= 1; ++q) rows.emplace_back(n, static_cast<std::uint8_t>(2 + q));
        } else if (k == n) {
            for (int p = 0; p <= 1; ++p) rows.emplace_back(n, static_cast<std::uint8_t>(p));
        } else {
            for (int p = 0; p <= 1; ++p) {
                for (int q = 0; q <= 1; ++q) {
                    std::vector<std::uint8_t> row(n);
                    for (int i = 0; i < k; ++i) row[i] = static_cast<std::uint8_t>(p);
                    for (int i = k; i < n; ++i) row[i] = static_cast<std::uint8_t>(2 + q);
                    rows.push_back(std::move(row));
                }
            }
        }
    }

    const int states = static_cast<int>(rows.size());
    std::vector<std::vector<int>> trans(states);
    for (int i = 0; i < states; ++i) {
        trans[i].reserve(states);
        for (int j = 0; j < states; ++j) {
            bool ok = true;
            for (int x = 0; x < n; ++x) {
                if (attacks(rows[i][x], rows[j][x], 0, 1)) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;
            for (int x = 0; x + 1 < n; ++x) {
                if (attacks(rows[i][x], rows[j][x + 1], 1, 1) ||
                    attacks(rows[i][x + 1], rows[j][x], -1, 1)) {
                    ok = false;
                    break;
                }
            }
            if (ok) trans[i].push_back(j);
        }
    }

    std::vector<u64> dp(states, 1), next(states, 0);
    for (int row = 1; row < n; ++row) {
        std::fill(next.begin(), next.end(), 0ULL);
        for (int i = 0; i < states; ++i) {
            if (dp[i] == 0) continue;
            for (int j : trans[i]) next[j] += dp[i];
        }
        dp.swap(next);
    }

    u64 total = 0;
    for (u64 v : dp) total += v;
    return total;
}

int C_mod(u64 n, const BlockFactorial& bf) {
    const int central = bf.central_binom(n);
    const int n_mod = static_cast<int>(n % static_cast<u64>(kMod));
    const int n2_mod = static_cast<int>((static_cast<u64>(n_mod) * n_mod) % kMod);
    const int correction = static_cast<int>(
        (3ULL * static_cast<u64>(n2_mod) + 2ULL * static_cast<u64>(n_mod) + 7ULL) % kMod);

    int ans = static_cast<int>((8ULL * static_cast<u64>(central)) % kMod);
    ans -= correction;
    if (ans < 0) ans += kMod;
    return ans;
}

std::vector<u64> fibonacci_upto_90() {
    std::vector<u64> f(kFibLimit + 1, 0);
    f[1] = 1;
    f[2] = 1;
    for (int i = 3; i <= kFibLimit; ++i) f[i] = f[i - 1] + f[i - 2];
    return f;
}

bool run_validations(const BlockFactorial& bf) {
    // Checkpoint 1: exact transfer counting on small n.
    struct ExactCheck {
        int n;
        u64 expected;
    };
    const std::vector<ExactCheck> exact_checks = {
        {1, 4ULL},     {2, 25ULL},     {3, 120ULL},    {4, 497ULL},     {5, 1924ULL},
        {6, 7265ULL},  {7, 27288ULL},  {8, 102745ULL}, {9, 388692ULL},  {10, 1477721ULL},
    };
    for (const auto& check : exact_checks) {
        const u64 got = exact_count_transfer(check.n);
        if (got != check.expected) {
            std::cerr << "Validation failed: exact C(" << check.n << ") = " << got
                      << ", expected " << check.expected << "\n";
            return false;
        }
        const int got_mod = C_mod(static_cast<u64>(check.n), bf);
        if (got_mod != static_cast<int>(check.expected % kMod)) {
            std::cerr << "Validation failed: formula C(" << check.n << ") mod p = " << got_mod
                      << ", expected " << (check.expected % kMod) << "\n";
            return false;
        }
    }

    // Checkpoint 2: binom_small for n<mod against Pascal DP.
    constexpr int kSmallN = 320;
    std::vector<int> row(kSmallN + 1, 0);
    row[0] = 1;
    for (int n = 0; n <= kSmallN; ++n) {
        if (n > 0) {
            for (int k = n; k >= 1; --k) {
                int v = row[k] + row[k - 1];
                if (v >= kMod) v -= kMod;
                row[k] = v;
            }
        }
        for (int k = 0; k <= n; ++k) {
            const int got = bf.binom_small(n, k);
            if (got != row[k]) {
                std::cerr << "Validation failed: binom_small(" << n << "," << k << ") = " << got
                          << ", expected " << row[k] << "\n";
                return false;
            }
        }
    }

    // Checkpoint 3: the statement examples.
    if (C_mod(1ULL, bf) != 4 || C_mod(2ULL, bf) != 25 || C_mod(10ULL, bf) != 1477721) {
        std::cerr << "Validation failed: sample values C(1), C(2), C(10) mismatch.\n";
        return false;
    }

    return true;
}

int solve_sum(const BlockFactorial& bf, unsigned threads) {
    const auto fib = fibonacci_upto_90();
    std::vector<int> values(kFibLimit + 1, 0);

    if (threads == 0) threads = 1;
    threads = std::min<unsigned>(threads, kFibLimit - 1);
    if (threads == 0) threads = 1;

    std::atomic<int> next_i(2);
    auto worker = [&]() {
        while (true) {
            const int i = next_i.fetch_add(1, std::memory_order_relaxed);
            if (i > kFibLimit) break;
            values[i] = C_mod(fib[i], bf);
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    int answer = 0;
    for (int i = 2; i <= kFibLimit; ++i) {
        answer += values[i];
        if (answer >= kMod) answer -= kMod;
    }
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;
    if (argc > 1) threads = static_cast<unsigned>(std::max(1, std::atoi(argv[1])));

    int block_size = kDefaultBlockSize;
    if (argc > 2) block_size = std::max(256, std::atoi(argv[2]));

    BlockFactorial bf(kMod, block_size);
    if (!run_validations(bf)) return 1;

    const int answer = solve_sum(bf, threads);
    std::cout << answer << "\n";
    return 0;
}

