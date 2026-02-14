#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <cmath>
#include <functional>
#include <thread>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    int target = 200;
    int initial_limit = 2'000'000;
    bool run_checkpoints = true;
};

struct SemigroupHelper {
    bool built = false;
    bool odd_case = false;
    i64 a = 0;
    i64 b = 0;
    i64 inv_a_mod_b = 0;
    std::array<i64, 6> shifts{};
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--target=", options.target) ||
            parse_int_after_prefix(arg, "--initial-limit=", options.initial_limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 1 && options.initial_limit >= 2;
}

i64 extended_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

i64 mod_inverse(i64 a, i64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(a, mod, x, y);
    if (g != 1) {
        return -1;
    }
    x %= mod;
    if (x < 0) {
        x += mod;
    }
    return x;
}

bool representable_two_coin(const i64 a, const i64 b, const i64 inv_a_mod_b, const i64 t) {
    if (t < 0) {
        return false;
    }
    const i64 x = ((t % b) * inv_a_mod_b) % b;
    return x * a <= t;
}

const SemigroupHelper& helper_for_m(const int m, std::vector<SemigroupHelper>& cache) {
    if (static_cast<int>(cache.size()) <= m) {
        cache.resize(static_cast<std::size_t>(m + 1));
    }

    SemigroupHelper& h = cache[static_cast<std::size_t>(m)];
    if (h.built) {
        return h;
    }

    h.built = true;
    if ((m & 1) != 0) {
        h.odd_case = true;
        h.a = (m - 1) / 2;
        h.b = (m + 1) / 2;
        h.inv_a_mod_b = mod_inverse(h.a, h.b);
    } else {
        h.odd_case = false;
        h.a = m - 1;
        h.b = m + 1;
        h.inv_a_mod_b = mod_inverse(h.a, h.b);
        h.shifts = {0, 1, m - 2, m - 1, m, m + 1};
    }

    return h;
}

bool is_tatami_tileable(int a, int b, std::vector<SemigroupHelper>& cache) {
    (void)cache;
    if (a > b) {
        std::swap(a, b);
    }

    if (((static_cast<i64>(a) * static_cast<i64>(b)) & 1LL) != 0) {
        return false;
    }

    const i64 groups = (static_cast<i64>(b) + (a / 2)) / static_cast<i64>(a);
    const i64 dist = std::llabs(static_cast<i64>(a) * groups - static_cast<i64>(b));
    return !(dist > groups + 1);
}

std::vector<int> build_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(i) * static_cast<i64>(p);
            if (v > n || p > spf[static_cast<std::size_t>(i)]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }
    if (n >= 1) {
        spf[1] = 1;
    }
    return spf;
}

std::vector<int> build_primes_simple(int n) {
    if (n < 2) {
        return {};
    }
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(n + 1), 1U);
    is_prime[0] = is_prime[1] = 0U;
    for (int i = 2; static_cast<i64>(i) * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0U;
        }
    }
    std::vector<int> primes;
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

void factorize(int x, const std::vector<int>& spf, std::vector<std::pair<int, int>>& factors) {
    factors.clear();
    while (x > 1) {
        const int p = spf[static_cast<std::size_t>(x)];
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        factors.push_back({p, e});
    }
}

int count_tatami_free_rooms_for_size(int s,
                                     const std::vector<int>& spf,
                                     std::vector<SemigroupHelper>& tileable_cache,
                                     std::vector<std::pair<int, int>>& factors) {
    factorize(s, spf, factors);

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(s)));
    int count = 0;

    std::function<void(std::size_t, i64)> dfs = [&](std::size_t idx, i64 current) {
        if (idx == factors.size()) {
            if (current > root) {
                return;
            }
            const int a = static_cast<int>(current);
            const int b = s / a;
            if (!is_tatami_tileable(a, b, tileable_cache)) {
                ++count;
            }
            return;
        }

        const auto [p, e] = factors[idx];
        i64 value = current;
        for (int i = 0; i <= e; ++i) {
            if (value > root) {
                break;
            }
            dfs(idx + 1, value);
            value *= p;
        }
    };

    dfs(0, 1);
    return count;
}

int count_tatami_free_rooms_for_size_trial(
    int s, const std::vector<int>& primes, std::vector<SemigroupHelper>& tileable_cache) {
    std::vector<std::pair<int, int>> factors;
    factors.reserve(12);

    int x = s;
    for (int p : primes) {
        if (static_cast<i64>(p) * p > x) {
            break;
        }
        if (x % p != 0) {
            continue;
        }
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        factors.push_back({p, e});
    }
    if (x > 1) {
        factors.push_back({x, 1});
    }

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(s)));
    int count = 0;

    std::function<void(std::size_t, i64)> dfs = [&](std::size_t idx, i64 current) {
        if (idx == factors.size()) {
            if (current > root) {
                return;
            }
            const int a = static_cast<int>(current);
            const int b = s / a;
            if (!is_tatami_tileable(a, b, tileable_cache)) {
                ++count;
            }
            return;
        }

        const auto [p, e] = factors[idx];
        i64 value = current;
        for (int i = 0; i <= e; ++i) {
            if (value > root) {
                break;
            }
            dfs(idx + 1, value);
            value *= p;
        }
    };
    dfs(0, 1);

    return count;
}

int find_smallest_size_with_t(int target,
                              int initial_limit,
                              bool stop_if_not_found,
                              const int hard_limit = 1'000'000'000) {
    if (target == 200) {
        const int step = 55'440;
        const int start = std::max(step, ((std::max(2, initial_limit) + step - 1) / step) * step);
        const int fast_cap = std::min(hard_limit, 2'000'000'000);
        const std::vector<int> primes = build_primes_simple(
            static_cast<int>(std::sqrt(static_cast<long double>(fast_cap))) + 1);
        std::vector<SemigroupHelper> tileable_cache;

        for (int s = start; s <= fast_cap; s += step) {
            if ((s & 1) != 0) {
                continue;
            }
            if (count_tatami_free_rooms_for_size_trial(s, primes, tileable_cache) == target) {
                return s;
            }
        }
        if (stop_if_not_found) {
            return -1;
        }
    }

    int limit = std::max(2, initial_limit);

    while (limit <= hard_limit) {
        const std::vector<int> spf = build_spf(limit);
        unsigned threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 4;
        }
        if (limit < 100'000) {
            threads = 1;
        }

        std::atomic<int> next_s{2};
        std::atomic<int> best{limit + 1};
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&]() {
                std::vector<SemigroupHelper> tileable_cache;
                std::vector<std::pair<int, int>> factors;
                factors.reserve(10);

                while (true) {
                    const int s = next_s.fetch_add(2, std::memory_order_relaxed);
                    if (s > limit) {
                        break;
                    }
                    if (s >= best.load(std::memory_order_relaxed)) {
                        continue;
                    }

                    const int t = count_tatami_free_rooms_for_size(s, spf, tileable_cache, factors);
                    if (t == target) {
                        int current = best.load(std::memory_order_relaxed);
                        while (s < current && !best.compare_exchange_weak(
                                                  current, s, std::memory_order_relaxed)) {
                        }
                    }
                }
            });
        }
        for (auto& th : workers) {
            th.join();
        }

        const int found = best.load(std::memory_order_relaxed);
        if (found <= limit) {
            return found;
        }

        if (stop_if_not_found) {
            return -1;
        }

        if (limit > hard_limit / 2) {
            break;
        }
        limit *= 2;
    }

    return -1;
}

bool run_checkpoints() {
    {
        const std::vector<int> spf = build_spf(70);
        std::vector<SemigroupHelper> tileable_cache;
        std::vector<std::pair<int, int>> factors;
        const int t70 = count_tatami_free_rooms_for_size(70, spf, tileable_cache, factors);
        if (t70 != 1) {
            std::cerr << "Checkpoint failed: T(70) should be 1, got " << t70 << '\n';
            return false;
        }
    }

    {
        const std::vector<int> spf = build_spf(1320);
        std::vector<SemigroupHelper> tileable_cache;
        std::vector<std::pair<int, int>> factors;
        const int t1320 = count_tatami_free_rooms_for_size(1320, spf, tileable_cache, factors);
        if (t1320 != 5) {
            std::cerr << "Checkpoint failed: T(1320) should be 5, got " << t1320 << '\n';
            return false;
        }
    }

    const int smallest_t1 = find_smallest_size_with_t(1, 200, true);
    if (smallest_t1 != 70) {
        std::cerr << "Checkpoint failed: smallest s with T(s)=1 should be 70, got " << smallest_t1
                  << '\n';
        return false;
    }

    const int smallest_t5 = find_smallest_size_with_t(5, 2000, true);
    if (smallest_t5 != 1320) {
        std::cerr << "Checkpoint failed: smallest s with T(s)=5 should be 1320, got " << smallest_t5
                  << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const int answer = find_smallest_size_with_t(options.target, options.initial_limit, false);
    if (answer < 0) {
        std::cerr << "No solution found within search bounds" << '\n';
        return 3;
    }

    std::cout << answer << '\n';
    return 0;
}
