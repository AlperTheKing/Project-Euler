#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }

    for (int i = 2; static_cast<long long>(i) * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

std::vector<std::vector<i64>> build_binom(int n) {
    std::vector<std::vector<i64>> c(n + 1, std::vector<i64>(n + 1, 0));
    for (int i = 0; i <= n; ++i) {
        c[i][0] = 1;
        c[i][i] = 1;
        for (int j = 1; j < i; ++j) {
            c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
        }
    }
    return c;
}

std::vector<std::vector<i64>> build_weights(int max_a, int max_b) {
    const int dim = std::max(max_a, max_b);
    const std::vector<std::vector<i64>> binom = build_binom(dim);

    std::vector<std::vector<i64>> w(max_a + 1, std::vector<i64>(max_b + 1, 0));

    for (int a = 0; a <= max_a; ++a) {
        for (int b = 0; b <= max_b; ++b) {
            i64 value = 0;
            for (int i = 0; i <= a; ++i) {
                for (int j = 0; j <= b; ++j) {
                    if ((i + 2 * j) % 3 != 0) {
                        continue;
                    }
                    i64 term = binom[a][i] * binom[b][j];
                    if ((i + j) & 1) {
                        term = -term;
                    }
                    value += term;
                }
            }
            if ((a + b) & 1) {
                value = -value;
            }
            w[a][b] = value;
        }
    }

    return w;
}

u64 solve(u64 N, int B) {
    std::vector<int> primes = primes_up_to(B);
    primes.erase(std::remove(primes.begin(), primes.end(), 3), primes.end());

    int max_a = 0;
    int max_b = 0;
    for (int p : primes) {
        if (p % 3 == 1) {
            ++max_a;
        } else {
            ++max_b;
        }
    }

    const std::vector<std::vector<i64>> w = build_weights(max_a, max_b);

    u64 answer = 0;

    auto dfs = [&](auto&& self, int idx, u64 product, int count_a, int count_b) -> void {
        if (idx == static_cast<int>(primes.size())) {
            const i64 weight = w[count_a][count_b];
            if (weight >= 0) {
                answer += static_cast<u64>(weight) * (N / product);
            } else {
                answer -= static_cast<u64>(-weight) * (N / product);
            }
            return;
        }

        self(self, idx + 1, product, count_a, count_b);

        const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(idx)]);
        if (product > N / p) {
            return;
        }

        if (primes[static_cast<std::size_t>(idx)] % 3 == 1) {
            self(self, idx + 1, product * p, count_a + 1, count_b);
        } else {
            self(self, idx + 1, product * p, count_a, count_b + 1);
        }
    };

    dfs(dfs, 0, 1, 0, 0);
    return answer;
}

bool run_checkpoints() {
    if (solve(10, 4) != 5ULL) {
        std::cerr << "Checkpoint failed: F(10,4) != 5" << '\n';
        return false;
    }
    if (solve(10, 10) != 3ULL) {
        std::cerr << "Checkpoint failed: F(10,10) != 3" << '\n';
        return false;
    }
    if (solve(100, 10) != 41ULL) {
        std::cerr << "Checkpoint failed: F(100,10) != 41" << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    std::cout << solve(1'000'000'000'000'000'000ULL, 120) << '\n';
    return 0;
}
