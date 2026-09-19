#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
constexpr u64 MOD = 1'234'567'891;
constexpr u64 PERIOD = MOD - 1;

constexpr bool is_prime(u64 value) {
    if (value < 2) return false;
    for (u64 d = 2; d * d <= value; ++d) {
        if (value % d == 0) return false;
    }
    return true;
}

static_assert(is_prime(MOD), "Exponent reduction requires a prime modulus");

u64 power(u64 base, u64 exponent) {
    u64 result = 1;
    while (exponent != 0) {
        if (exponent & 1) result = result * base % MOD;
        base = base * base % MOD;
        exponent >>= 1;
    }
    return result;
}

u64 add(u64 a, u64 b) {
    const u64 sum = a + b;
    return sum >= PERIOD ? sum - PERIOD : sum;
}

struct Counts {
    u64 partitions;
    std::vector<u64> hooks;
};

Counts count_hooks(int m, int n) {
    if (m < 1 || n < 0) throw std::invalid_argument("Invalid penguin or step count");
    if (n == 0) return {1, {0}};
    m = std::min(m, n);
    const int baseline = m * (m - 1) / 2;
    const int degree = n + baseline;
    std::vector<u64> partitions(degree + 1), quotient(degree + 1);
    Counts result{0, std::vector<u64>(n + 1)};
    partitions[0] = 1;

    // P_r(q) = product_(i=1..r) 1/(1-q^i); beta numbers encode the hooks.
    for (int r = 0; r < m; ++r) {
        const int j = m - r;
        const int shift = baseline - r * (r - 1) / 2;
        for (int s = 0; s <= n + shift; ++s) {
            quotient[s] = add(partitions[s], s >= j ? quotient[s - j] : 0);
        }

        // H_h = sum_j (-1)^(j-1) sum_(t=1..j) [q^(n+shift-t*h)] P_r/(1-q^j).
        for (int h = 1; h <= n; ++h) {
            u64 contribution = 0;
            for (int t = 1, s = n + shift - h; t <= j && s >= 0; ++t, s -= h) {
                contribution += quotient[s];
            }
            contribution %= PERIOD;
            result.hooks[h] = add(result.hooks[h],
                j % 2 != 0 ? contribution : (PERIOD - contribution) % PERIOD);
        }

        const int part = r + 1;
        for (int s = part; s <= degree; ++s) {
            partitions[s] = add(partitions[s], partitions[s - part]);
        }
    }
    result.partitions = partitions[n];
    return result;
}

u64 product_from_hooks(const Counts& counts) {
    u64 factorial = 1, denominator = 1;
    for (std::size_t h = 1; h < counts.hooks.size(); ++h) {
        factorial = factorial * h % MOD;
        denominator = denominator * power(h, counts.hooks[h]) % MOD;
    }
    // Each endpoint contributes n! divided by the product of its hook lengths.
    return power(factorial, counts.partitions) * power(denominator, MOD - 2) % MOD;
}

u64 solve(int m, int n) {
    return product_from_hooks(count_hooks(m, n));
}

void require(bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

void check_walks(int m, int max_steps) {
    std::vector<int> initial(m);
    std::iota(initial.begin(), initial.end(), 1);
    std::map<std::vector<int>, u64> states{{initial, 1}};
    for (int n = 0; n <= max_steps; ++n) {
        u64 expected = 1;
        std::vector<u64> hooks(n + 1);
        for (const auto& [positions, ways] : states) {
            expected = expected * (ways % MOD) % MOD;
            std::vector<int> shape(m);
            for (int i = 0; i < m; ++i) shape[i] = positions[m - 1 - i] - (m - i);
            for (int row = 0; row < m; ++row) {
                for (int col = 0; col < shape[row]; ++col) {
                    int hook = shape[row] - col;
                    for (int below = row + 1; below < m; ++below) {
                        if (shape[below] > col) ++hook;
                    }
                    ++hooks[hook];
                }
            }
        }
        const Counts actual = count_hooks(m, n);
        const std::string label = "m=" + std::to_string(m) + ", n=" + std::to_string(n);
        require(actual.partitions == states.size(), "endpoint count for " + label);
        require(actual.hooks == hooks, "hook counts for " + label);
        require(product_from_hooks(actual) == expected, "legal-walk product for " + label);
        if (n == max_steps) break;
        std::map<std::vector<int>, u64> next;
        for (const auto& [positions, ways] : states) {
            for (int i = 0; i < m; ++i) {
                if (i + 1 < m && positions[i] + 1 == positions[i + 1]) continue;
                auto moved = positions;
                ++moved[i];
                next[moved] += ways;
            }
        }
        states = std::move(next);
    }
}

void run_tests() {
    require(solve(2, 4) == 6, "F(2,4)");
    require(solve(3, 6) == 180'000, "F(3,6)");
    require(solve(5, 10) == 411'456'133, "F(5,10)");
    for (int m = 1; m <= 16; ++m) check_walks(m, 16);
    const Counts target = count_hooks(150, 300);
    u64 total = 0;
    for (const u64 count : target.hooks) total = add(total, count);
    require(total == 300 * target.partitions % PERIOD, "total target hook count");
    std::cout << "All checks passed.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
            return EXIT_SUCCESS;
        }
        if (argc != 1) throw std::invalid_argument("Usage: Euler1010 [--self-test]");
        std::cout << solve(150, 300) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
