#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

namespace {

long double combination_ld(int n, int k) {
    if (k < 0 || k > n) return 0.0L;
    k = std::min(k, n - k);
    long double value = 1.0L;
    for (int i = 1; i <= k; ++i) {
        value *= static_cast<long double>(n - k + i);
        value /= static_cast<long double>(i);
    }
    return value;
}

long double probability_exact_marked_fixed(int n, int marked, int fixed_marked) {
    if (n < 0 || marked < 0 || marked > n) return 0.0L;
    if (fixed_marked < 0 || fixed_marked > marked) return 0.0L;

    const int bad_marked = marked - fixed_marked;

    long double term = 1.0L;
    for (int i = 0; i < fixed_marked; ++i) {
        term /= static_cast<long double>(n - i);
    }

    long double alternating_sum = term;
    for (int j = 0; j < bad_marked; ++j) {
        const long double numer = -static_cast<long double>(bad_marked - j);
        const long double denom = static_cast<long double>(j + 1) *
                                  static_cast<long double>(n - fixed_marked - j);
        term *= numer / denom;
        alternating_sum += term;
    }

    return combination_ld(marked, fixed_marked) * alternating_sum;
}

__int128 combination_i128(int n, int k) {
    if (k < 0 || k > n) return 0;
    k = std::min(k, n - k);
    __int128 value = 1;
    for (int i = 1; i <= k; ++i) {
        value = value * static_cast<__int128>(n - k + i) / static_cast<__int128>(i);
    }
    return value;
}

__int128 factorial_i128(int n) {
    __int128 value = 1;
    for (int i = 2; i <= n; ++i) {
        value *= static_cast<__int128>(i);
    }
    return value;
}

__int128 count_exact_marked_fixed_inclusion(int n, int marked, int fixed_marked) {
    if (n < 0 || marked < 0 || marked > n) return 0;
    if (fixed_marked < 0 || fixed_marked > marked) return 0;

    const int bad_marked = marked - fixed_marked;
    __int128 alternating_sum = 0;
    for (int j = 0; j <= bad_marked; ++j) {
        __int128 term = combination_i128(bad_marked, j) * factorial_i128(n - fixed_marked - j);
        if ((j & 1) == 0) {
            alternating_sum += term;
        } else {
            alternating_sum -= term;
        }
    }
    return combination_i128(marked, fixed_marked) * alternating_sum;
}

std::uint64_t brute_count_exact_marked_fixed(int n, int marked, int fixed_marked, int threads) {
    if (threads < 1) threads = 1;
    if (n <= 0) return (marked == 0 && fixed_marked == 0) ? 1ULL : 0ULL;

    std::atomic<int> next_first{0};
    std::vector<std::uint64_t> partial(static_cast<std::size_t>(threads), 0ULL);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            std::uint64_t local = 0;
            while (true) {
                const int first = next_first.fetch_add(1, std::memory_order_relaxed);
                if (first >= n) break;

                std::vector<int> rest;
                rest.reserve(static_cast<std::size_t>(n - 1));
                for (int v = 0; v < n; ++v) {
                    if (v != first) rest.push_back(v);
                }

                do {
                    int fixed = 0;
                    if (marked > 0 && first == 0) ++fixed;
                    for (int pos = 1; pos < marked; ++pos) {
                        if (rest[static_cast<std::size_t>(pos - 1)] == pos) ++fixed;
                    }
                    if (fixed == fixed_marked) ++local;
                } while (std::next_permutation(rest.begin(), rest.end()));
            }
            partial[static_cast<std::size_t>(t)] = local;
        });
    }

    for (std::thread& th : pool) {
        th.join();
    }

    std::uint64_t total = 0;
    for (std::uint64_t value : partial) total += value;
    return total;
}

bool validate() {
    struct TinyCase {
        int n;
        int marked;
        int fixed_marked;
    };

    const std::vector<TinyCase> tiny_cases = {
        {7, 3, 1},
        {8, 4, 2},
        {9, 5, 1},
        {10, 4, 2},
    };

    for (const TinyCase& c : tiny_cases) {
        const __int128 exact_count = count_exact_marked_fixed_inclusion(c.n, c.marked, c.fixed_marked);
        const std::uint64_t brute_count = brute_count_exact_marked_fixed(c.n, c.marked, c.fixed_marked, 1);
        if (exact_count != static_cast<__int128>(brute_count)) {
            std::cerr << "Validation failed (inclusion vs brute) at n=" << c.n
                      << ", marked=" << c.marked << ", fixed=" << c.fixed_marked
                      << ": inclusion=" << static_cast<long long>(exact_count)
                      << ", brute=" << brute_count << "\n";
            return false;
        }

        const long double prob_formula = probability_exact_marked_fixed(c.n, c.marked, c.fixed_marked);
        const long double prob_exact = static_cast<long double>(brute_count) /
                                       static_cast<long double>(factorial_i128(c.n));
        if (std::fabsl(prob_formula - prob_exact) > 1e-16L) {
            std::cerr << "Validation failed (probability mismatch) at n=" << c.n
                      << ", marked=" << c.marked << ", fixed=" << c.fixed_marked
                      << ": formula=" << std::setprecision(22) << prob_formula
                      << ", exact=" << prob_exact << "\n";
            return false;
        }
    }

    const int n_mass = 20;
    const int marked_mass = 7;
    long double mass = 0.0L;
    for (int fixed = 0; fixed <= marked_mass; ++fixed) {
        mass += probability_exact_marked_fixed(n_mass, marked_mass, fixed);
    }
    if (std::fabsl(mass - 1.0L) > 1e-16L) {
        std::cerr << "Validation failed (total mass): got=" << std::setprecision(22) << mass << "\n";
        return false;
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 2;
    const int threads = static_cast<int>(std::min<unsigned>(hw, 8));
    if (threads > 1) {
        const TinyCase thread_case{10, 4, 2};
        const std::uint64_t single =
            brute_count_exact_marked_fixed(thread_case.n, thread_case.marked, thread_case.fixed_marked, 1);
        const std::uint64_t multi =
            brute_count_exact_marked_fixed(thread_case.n, thread_case.marked, thread_case.fixed_marked, threads);
        if (single != multi) {
            std::cerr << "Validation failed (thread consistency): single=" << single
                      << ", multi=" << multi << "\n";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!validate()) return 1;

    int n = 100;
    int marked = 25;
    int misplaced_marked = 22;

    if (argc > 1) n = std::max(1, std::atoi(argv[1]));
    if (argc > 2) marked = std::clamp(std::atoi(argv[2]), 0, n);
    if (argc > 3) misplaced_marked = std::clamp(std::atoi(argv[3]), 0, marked);

    const int fixed_marked = marked - misplaced_marked;
    const long double answer = probability_exact_marked_fixed(n, marked, fixed_marked);
    std::cout << std::fixed << std::setprecision(12) << answer << '\n';
    return 0;
}
