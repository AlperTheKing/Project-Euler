#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;

static constexpr i64 MOD = 1'020'340'567LL;

static i64 solve(int N) {
    std::vector<int> lp(N + 1, 0);
    std::vector<int> primes;
    primes.reserve(N / 10);

    std::vector<int> mu(N + 1, 0);
    mu[1] = 1;

    for (int i = 2; i <= N; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            const i64 v = 1LL * i * p;
            if (v > N) {
                break;
            }
            lp[static_cast<int>(v)] = p;
            if (p == lp[i]) {
                mu[static_cast<int>(v)] = 0;
                break;
            }
            mu[static_cast<int>(v)] = -mu[i];
        }
    }

    std::vector<int> pref_mu(N + 1, 0);
    for (int i = 1; i <= N; ++i) {
        pref_mu[i] = pref_mu[i - 1] + mu[i];
    }

    std::vector<i64> pow2(N + 1, 1);
    for (int i = 1; i <= N; ++i) {
        pow2[i] = (pow2[i - 1] * 2LL) % MOD;
    }

    i64 ans = 0;
    int i = 1;
    while (i <= N) {
        const int q = N / i;
        const int j = N / q;

        const i64 block_mu = static_cast<i64>(pref_mu[j] - pref_mu[i - 1]);
        i64 term = block_mu % MOD;
        if (term < 0) {
            term += MOD;
        }

        ans += (term * pow2[q]) % MOD;
        ans %= MOD;

        i = j + 1;
    }

    if (ans < 0) {
        ans += MOD;
    }
    return ans;
}

int main() {
    assert(solve(1) == 2);
    assert(solve(2) == 2);
    assert(solve(3) == 4);

    std::cout << solve(10'000'000) << '\n';
    return 0;
}
