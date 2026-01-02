#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

using int64 = long long;

constexpr int64 MOD = 1000000123LL;

static int64 mod_pow(int64 base, int64 exp) {
    int64 result = 1 % MOD;
    base %= MOD;
    while (exp > 0) {
        if (exp & 1LL) result = (result * base) % MOD;
        base = (base * base) % MOD;
        exp >>= 1LL;
    }
    return result;
}

static void build_factorials(int n, std::vector<int64>& fact, std::vector<int64>& inv_fact) {
    fact.assign(n + 1, 0);
    inv_fact.assign(n + 1, 0);
    fact[0] = 1;
    for (int i = 1; i <= n; ++i) {
        fact[i] = (fact[i - 1] * i) % MOD;
    }
    inv_fact[n] = mod_pow(fact[n], MOD - 2);
    for (int i = n; i >= 1; --i) {
        inv_fact[i - 1] = (inv_fact[i] * i) % MOD;
    }
}

static std::vector<int64> build_pow_invfact(const std::vector<int64>& inv_fact, int n, int r) {
    std::vector<int64> pow_invfact(n + 1, 0);
    for (int i = 0; i <= n; ++i) {
        pow_invfact[i] = mod_pow(inv_fact[i], r);
    }
    return pow_invfact;
}

static int64 compute_core_with_workspace(int n, int k, const int64* pow_invfact, std::vector<int64>& dp) {
    const int blocks = n / k;
    const int rem = n - blocks * k;

    dp[0] = 1;
    for (int i = 1; i <= blocks; ++i) {
        int64 sum = 0;
        int sign = 1;
        int idx = k;
        for (int j = 1; j <= i; ++j) {
            int64 term = (dp[i - j] * pow_invfact[idx]) % MOD;
            sum += sign * term;
            sign = -sign;
            idx += k;
        }
        sum %= MOD;
        if (sum < 0) sum += MOD;
        dp[i] = sum;
    }

    if (rem == 0) {
        return dp[blocks];
    }

    int64 total = 0;
    int sign = 1;
    int idx = rem;
    for (int j = 0; j <= blocks; ++j) {
        int64 term = (dp[blocks - j] * pow_invfact[idx]) % MOD;
        total += sign * term;
        sign = -sign;
        idx += k;
    }
    total %= MOD;
    if (total < 0) total += MOD;
    return total;
}

static int64 compute_P(int n, int r, int k, const std::vector<int64>& fact,
                       const std::vector<int64>& inv_fact) {
    std::vector<int64> pow_invfact = build_pow_invfact(inv_fact, n, r);
    std::vector<int64> dp(n + 1, 0);
    int64 core = compute_core_with_workspace(n, k, pow_invfact.data(), dp);
    int64 factor = mod_pow(fact[n], r);
    return (core * factor) % MOD;
}

static int64 compute_Q(int n, const std::vector<int64>& fact, const std::vector<int64>& inv_fact,
                       unsigned threads) {
    std::vector<int64> pow_invfact = build_pow_invfact(inv_fact, n, n);
    const int64* w = pow_invfact.data();

    if (threads == 0) threads = 1;
    threads = std::min<unsigned>(threads, static_cast<unsigned>(n));

    std::atomic<int> next_k(1);
    std::vector<int64> partial(threads, 0);

    auto worker = [&](unsigned tid) {
        std::vector<int64> dp(n + 1, 0);
        int64 local = 0;
        while (true) {
            int k = next_k.fetch_add(1, std::memory_order_relaxed);
            if (k > n) break;

            int64 core = compute_core_with_workspace(n, k, w, dp);
            local += core;
        }
        partial[tid] = local % MOD;
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for (auto& th : pool) th.join();

    int64 total = 0;
    for (int64 v : partial) {
        total += v;
    }
    total %= MOD;
    int64 factor = mod_pow(fact[n], n);
    return (total * factor) % MOD;
}

static bool run_validations(const std::vector<int64>& fact, const std::vector<int64>& inv_fact) {
    struct CheckP {
        int k;
        int r;
        int n;
        int64 expected;
    };
    const std::vector<CheckP> checks = {
        {1, 2, 3, 19},
        {2, 4, 6, 65508751},
        {7, 5, 30, 161858102},
    };

    for (const auto& check : checks) {
        int64 got = compute_P(check.n, check.r, check.k, fact, inv_fact);
        if (got != check.expected) {
            std::cerr << "Validation failed for P(" << check.k << ", " << check.r << ", "
                      << check.n << "): got=" << got << " expected=" << check.expected << "\n";
            return false;
        }
    }

    const int64 expected_q5 = 879391168;
    int64 q5 = compute_Q(5, fact, inv_fact, 1);
    if (q5 != expected_q5) {
        std::cerr << "Validation failed for Q(5): got=" << q5 << " expected=" << expected_q5
                  << "\n";
        return false;
    }

    const int64 expected_q50 = 819573537;
    int64 q50 = compute_Q(50, fact, inv_fact, 1);
    if (q50 != expected_q50) {
        std::cerr << "Validation failed for Q(50): got=" << q50 << " expected=" << expected_q50
                  << "\n";
        return false;
    }

    return true;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n = 50000;
    if (argc > 1) n = std::atoi(argv[1]);
    if (n <= 0) {
        std::cerr << "n must be positive.\n";
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;
    if (argc > 2) threads = static_cast<unsigned>(std::max(1, std::atoi(argv[2])));

    int max_n = std::max(n, 50);
    std::vector<int64> fact;
    std::vector<int64> inv_fact;
    build_factorials(max_n, fact, inv_fact);

    if (!run_validations(fact, inv_fact)) return 1;

    int64 answer = compute_Q(n, fact, inv_fact, threads);
    std::cout << answer << "\n";
    return 0;
}
