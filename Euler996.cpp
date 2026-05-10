#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr int TARGET_N = 123;
constexpr u64 TARGET_K = 4'567'891ULL;
constexpr u64 TARGET_E = TARGET_K / 2;
constexpr u64 MOD = 1'234'567'891ULL;

std::vector<std::vector<u64>> binomial_table(const int n, const int r, const u64 mod) {
    std::vector<std::vector<u64>> c(static_cast<std::size_t>(n + 1),
                                    std::vector<u64>(static_cast<std::size_t>(r + 1), 0));
    c[0][0] = 1;
    for (int i = 1; i <= n; ++i) {
        c[static_cast<std::size_t>(i)][0] = 1;
        const int hi = std::min(i, r);
        for (int j = 1; j <= hi; ++j) {
            const u64 value = c[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)] +
                              c[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)];
            c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                mod == 0 ? value : value % mod;
        }
    }
    return c;
}

u64 choose_from_table(const std::vector<std::vector<u64>>& c, const int n, const int r) {
    if (n < 0 || r < 0 || n < r) {
        return 0;
    }
    if (n >= static_cast<int>(c.size()) || r >= static_cast<int>(c[0].size())) {
        return 0;
    }
    return c[static_cast<std::size_t>(n)][static_cast<std::size_t>(r)];
}

std::vector<std::vector<u64>> run_count_table(const int n, const int emax, const u64 mod) {
    const auto comb = binomial_table(2 * emax + 1, n, mod);
    std::vector<std::vector<u64>> runs(static_cast<std::size_t>(n + 1),
                                       std::vector<u64>(static_cast<std::size_t>(emax + 1), 0));

    for (int length = 2; length <= n; ++length) {
        for (int e = 1; e <= emax; ++e) {
            const u64 all_positive = choose_from_table(comb, 2 * e - 1, length - 1);
            const u64 too_large = choose_from_table(comb, e - 1, length - 1);
            if (mod == 0) {
                runs[static_cast<std::size_t>(length)][static_cast<std::size_t>(e)] =
                    all_positive - static_cast<u64>(length) * too_large;
            } else {
                const u64 bad = static_cast<u64>(length) * too_large % mod;
                runs[static_cast<std::size_t>(length)][static_cast<std::size_t>(e)] =
                    (all_positive + mod - bad) % mod;
            }
        }
    }

    return runs;
}

std::vector<u64> cumulative_counts(const int n, const int emax, const u64 mod) {
    const auto runs = run_count_table(n, emax, mod);
    std::vector<std::vector<u64>> total(static_cast<std::size_t>(n + 1),
                                        std::vector<u64>(static_cast<std::size_t>(emax + 1), 0));
    std::vector<std::vector<u64>> zero(static_cast<std::size_t>(n + 1),
                                       std::vector<u64>(static_cast<std::size_t>(emax + 1), 0));
    total[0][0] = 1;
    zero[0][0] = 1;

    for (int i = 1; i <= n; ++i) {
        zero[static_cast<std::size_t>(i)] = total[static_cast<std::size_t>(i - 1)];
        std::vector<u64> row = zero[static_cast<std::size_t>(i)];

        for (int length = 2; length <= i; ++length) {
            const auto& source = zero[static_cast<std::size_t>(i - length)];
            const auto& run = runs[static_cast<std::size_t>(length)];
            for (int before = 0; before <= emax; ++before) {
                const u64 source_value = source[static_cast<std::size_t>(before)];
                if (source_value == 0) {
                    continue;
                }
                for (int used = 1; before + used <= emax; ++used) {
                    const u64 run_value = run[static_cast<std::size_t>(used)];
                    if (run_value == 0) {
                        continue;
                    }
                    u64& cell = row[static_cast<std::size_t>(before + used)];
                    if (mod == 0) {
                        cell += source_value * run_value;
                    } else {
                        cell = (cell + source_value * run_value) % mod;
                    }
                }
            }
        }

        total[static_cast<std::size_t>(i)] = std::move(row);
    }

    std::vector<u64> cumulative(static_cast<std::size_t>(emax + 1), 0);
    u64 sum = 0;
    for (int e = 0; e <= emax; ++e) {
        if (mod == 0) {
            sum += total[static_cast<std::size_t>(n)][static_cast<std::size_t>(e)];
        } else {
            sum = (sum + total[static_cast<std::size_t>(n)][static_cast<std::size_t>(e)]) % mod;
        }
        cumulative[static_cast<std::size_t>(e)] = sum;
    }
    return cumulative;
}

u64 binomial_large_mod(const u64 n, int r) {
    if (r < 0 || n < static_cast<u64>(r)) {
        return 0;
    }
    if (static_cast<u64>(r) > n - static_cast<u64>(r)) {
        r = static_cast<int>(n - static_cast<u64>(r));
    }

    std::vector<u64> numerator(static_cast<std::size_t>(r));
    for (int i = 0; i < r; ++i) {
        numerator[static_cast<std::size_t>(i)] = n - static_cast<u64>(r) + 1 + static_cast<u64>(i);
    }

    for (int d = 2; d <= r; ++d) {
        u64 remaining = static_cast<u64>(d);
        for (u64& value : numerator) {
            const u64 g = std::gcd(value, remaining);
            if (g == 1) {
                continue;
            }
            value /= g;
            remaining /= g;
            if (remaining == 1) {
                break;
            }
        }
        assert(remaining == 1);
    }

    u64 result = 1;
    for (const u64 value : numerator) {
        result = result * (value % MOD) % MOD;
    }
    return result;
}

u64 polynomial_tail_value(const int n, const u64 e) {
    const int start = n - 1;
    const int degree = n;
    const int emax = start + degree;
    const std::vector<u64> cumulative = cumulative_counts(n, emax, MOD);

    std::vector<u64> current;
    current.reserve(static_cast<std::size_t>(degree + 1));
    for (int i = 0; i <= degree; ++i) {
        current.push_back(cumulative[static_cast<std::size_t>(start + i)]);
    }

    std::vector<u64> differences;
    differences.reserve(static_cast<std::size_t>(degree + 1));
    for (int r = 0; r <= degree; ++r) {
        differences.push_back(current[0]);
        std::vector<u64> next;
        next.reserve(current.size() - 1);
        for (std::size_t i = 0; i + 1 < current.size(); ++i) {
            next.push_back((current[i + 1] + MOD - current[i]) % MOD);
        }
        current = std::move(next);
    }

    const u64 x = e - static_cast<u64>(start);
    u64 result = 0;
    for (int r = 0; r <= degree; ++r) {
        result = (result + differences[static_cast<std::size_t>(r)] * binomial_large_mod(x, r)) % MOD;
    }
    return result;
}

u64 brute_count(const int n, const int emax) {
    using State = std::pair<std::vector<int>, std::vector<int>>;

    std::vector<int> identity(static_cast<std::size_t>(n));
    std::iota(identity.begin(), identity.end(), 0);

    std::set<State> states;
    states.insert({identity, std::vector<int>(static_cast<std::size_t>(n), 0)});

    std::set<std::vector<int>> results;
    results.insert(std::vector<int>(static_cast<std::size_t>(n), 0));

    for (int step = 1; step <= 2 * emax; ++step) {
        std::set<State> next;
        for (const auto& [perm, count] : states) {
            for (int p = 0; p + 1 < n; ++p) {
                std::vector<int> next_perm = perm;
                std::vector<int> next_count = count;
                ++next_count[static_cast<std::size_t>(next_perm[static_cast<std::size_t>(p + 1)])];
                std::swap(next_perm[static_cast<std::size_t>(p)], next_perm[static_cast<std::size_t>(p + 1)]);
                if (next_perm == identity) {
                    results.insert(next_count);
                }
                next.insert({std::move(next_perm), std::move(next_count)});
            }
        }
        states = std::move(next);
    }

    return static_cast<u64>(results.size());
}

void run_checkpoints() {
    for (int n = 2; n <= 5; ++n) {
        const std::vector<u64> values = cumulative_counts(n, 4, 0);
        for (int e = 0; e <= 4; ++e) {
            assert(values[static_cast<std::size_t>(e)] == brute_count(n, e));
        }
    }

    assert(cumulative_counts(3, 2, 0)[2] == 8);
    assert(cumulative_counts(12, 17, 0)[17] == 2'457'178'250ULL);

    for (int n = 2; n <= 8; ++n) {
        const int start = n - 1;
        const int degree = n;
        const int emax = start + degree + 5;
        const std::vector<u64> values = cumulative_counts(n, emax, 0);

        std::vector<u64> current;
        for (int i = 0; i <= degree; ++i) {
            current.push_back(values[static_cast<std::size_t>(start + i)]);
        }

        std::vector<u64> differences;
        for (int r = 0; r <= degree; ++r) {
            differences.push_back(current[0]);
            std::vector<u64> next;
            for (std::size_t i = 0; i + 1 < current.size(); ++i) {
                next.push_back(current[i + 1] - current[i]);
            }
            current = std::move(next);
        }

        for (int e = start; e <= emax; ++e) {
            u64 predicted = 0;
            u64 choose = 1;
            for (int r = 0; r <= degree; ++r) {
                if (r > 0) {
                    choose = choose * static_cast<u64>(e - start - r + 1) / static_cast<u64>(r);
                }
                predicted += differences[static_cast<std::size_t>(r)] * choose;
            }
            assert(predicted == values[static_cast<std::size_t>(e)]);
        }
    }
}

}  // namespace

int main() {
    run_checkpoints();
    std::cout << polynomial_tail_value(TARGET_N, TARGET_E) << '\n';
    return 0;
}
