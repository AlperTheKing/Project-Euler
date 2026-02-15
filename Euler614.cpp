#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

using i64 = long long;

static constexpr int MOD = 1'000'000'007;

static std::vector<int> distinct_parts(int limit) {
    std::vector<int> q(static_cast<std::size_t>(limit + 1), 0);
    q[0] = 1;

    std::vector<int> pent_sign(static_cast<std::size_t>(limit + 1), 0);
    const int lim = static_cast<int>(std::sqrt(static_cast<double>(limit)));
    for (int j = 0; j <= lim; ++j) {
        const int sign = (j & 1) ? -1 : 1;
        const i64 g1 = static_cast<i64>(j) * (3LL * j + 1) / 2;
        const i64 g2 = static_cast<i64>(j) * (3LL * j - 1) / 2;
        if (g1 <= limit) pent_sign[static_cast<std::size_t>(g1)] = sign;
        if (g2 <= limit) pent_sign[static_cast<std::size_t>(g2)] = sign;
    }

    std::vector<int> squares;
    squares.reserve(static_cast<std::size_t>(lim));
    for (int k = 1; k <= lim; ++k) squares.push_back(k * k);

    for (int n = 1; n <= limit; ++n) {
        int val = pent_sign[static_cast<std::size_t>(n)];
        int sign = 1;
        for (int k = 0; k < static_cast<int>(squares.size()); ++k) {
            const int sq = squares[static_cast<std::size_t>(k)];
            if (sq > n) break;
            int term = q[static_cast<std::size_t>(n - sq)] << 1;
            if (term >= MOD) term -= MOD;
            if (sign > 0) {
                val += term;
                if (val >= MOD) val -= MOD;
            } else {
                val -= term;
                if (val < 0) val += MOD;
            }
            sign = -sign;
        }
        q[static_cast<std::size_t>(n)] = val;
    }
    return q;
}

static std::vector<int> self_conjugate(int limit) {
    const int kmax = static_cast<int>(std::sqrt(static_cast<double>(limit)));
    std::vector<int> q_odd(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> p(static_cast<std::size_t>(limit / 2 + 1), 0);
    q_odd[0] = 1;
    p[0] = 1;

    int* qptr = q_odd.data();
    int* pptr = p.data();

    for (int k = 1; k <= kmax; ++k) {
        const int k_sq = k * k;
        const int rem = limit - k_sq;
        const int vmax_small = std::min(k - 1, rem / 2);
        int idx = k_sq;
        for (int v = 0; v <= vmax_small; ++v, idx += 2) {
            int x = qptr[idx] + pptr[v];
            if (x >= MOD) x -= MOD;
            qptr[idx] = x;
        }

        const int vmax = rem / 2;
        int idx2 = k_sq + 2 * k;
        int* p_cur = pptr + k;
        int* p_prev = pptr;
        for (int v = k; v <= vmax; ++v, idx2 += 2, ++p_cur, ++p_prev) {
            int pv = *p_cur + *p_prev;
            if (pv >= MOD) pv -= MOD;
            *p_cur = pv;
            int x = qptr[idx2] + pv;
            if (x >= MOD) x -= MOD;
            qptr[idx2] = x;
        }
    }

    return q_odd;
}

static int solve614(int limit) {
    struct Task {
        int limit;
        std::vector<int>* out;
    };

    auto worker = [](void* arg) -> void* {
        auto* task = static_cast<Task*>(arg);
        *task->out = distinct_parts(task->limit);
        return nullptr;
    };

    std::vector<int> q;
    Task task{limit / 4, &q};
    pthread_t thread{};
    pthread_create(&thread, nullptr, worker, &task);

    std::vector<int> q_odd = self_conjugate(limit);
    pthread_join(thread, nullptr);

    std::vector<int> s_odd = q_odd;
    for (int i = 1; i <= limit; ++i) {
        int x = s_odd[static_cast<std::size_t>(i)] + s_odd[static_cast<std::size_t>(i - 1)];
        if (x >= MOD) x -= MOD;
        s_odd[static_cast<std::size_t>(i)] = x;
    }

    int res = MOD - 1;
    for (int i = 0; i <= limit; i += 4) {
        res = static_cast<int>((res + 1LL * s_odd[static_cast<std::size_t>(limit - i)] *
                                           q[static_cast<std::size_t>(i / 4)]) %
                               MOD);
    }
    return res;
}

static int brute_cumulative(int limit) {
    std::vector<int> dp(static_cast<std::size_t>(limit + 1), 0);
    dp[0] = 1;
    for (int p = 1; p <= limit; ++p) {
        if ((p & 3) == 2) continue;
        for (int s = limit; s >= p; --s) {
            int x = dp[static_cast<std::size_t>(s)] + dp[static_cast<std::size_t>(s - p)];
            if (x >= MOD) x -= MOD;
            dp[static_cast<std::size_t>(s)] = x;
        }
    }
    int sum = 0;
    for (int s = 1; s <= limit; ++s) {
        sum += dp[static_cast<std::size_t>(s)];
        if (sum >= MOD) sum -= MOD;
    }
    return sum;
}

int main(int argc, char** argv) {
    int n = 10'000'000;
    if (argc >= 2) n = std::stoi(argv[1]);

    assert(solve614(100) == brute_cumulative(100));
    assert(solve614(500) == brute_cumulative(500));
    assert(solve614(1000) == brute_cumulative(1000));

    std::cout << solve614(n) << '\n';
    return 0;
}
