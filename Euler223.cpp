#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

namespace {

constexpr int kDefaultLimit = 25'000'000;

struct Factor {
    int p;
    int e;
};

std::vector<int> build_spf(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            long long v = 1LL * p * i;
            if (v > n || p > spf[i]) break;
            spf[static_cast<std::size_t>(v)] = p;
        }
    }
    if (n >= 1) spf[1] = 1;
    return spf;
}

void factorize_int(int x, const std::vector<int>& spf, std::vector<Factor>& out) {
    out.clear();
    while (x > 1) {
        int p = spf[static_cast<std::size_t>(x)];
        int e = 0;
        do {
            x /= p;
            ++e;
        } while (x > 1 && spf[static_cast<std::size_t>(x)] == p);
        out.push_back({p, e});
    }
}

void merge_factors(const std::vector<Factor>& a,
                   const std::vector<Factor>& b,
                   std::vector<Factor>& out) {
    out.clear();
    out.reserve(a.size() + b.size());

    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a.size() || j < b.size()) {
        if (j == b.size() || (i < a.size() && a[i].p < b[j].p)) {
            out.push_back(a[i++]);
        } else if (i == a.size() || b[j].p < a[i].p) {
            out.push_back(b[j++]);
        } else {
            out.push_back({a[i].p, a[i].e + b[j].e});
            ++i;
            ++j;
        }
    }
}

std::uint64_t count_for_a(int a,
                          int perimeter_limit,
                          const std::vector<int>& spf,
                          std::vector<Factor>& left_factors,
                          std::vector<Factor>& right_factors,
                          std::vector<Factor>& merged_factors,
                          std::vector<std::uint64_t>& divisors) {
    const std::uint64_t n = static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(a) - 1ULL;
    if (n == 0) return 0;

    const std::uint64_t denom = static_cast<std::uint64_t>(perimeter_limit - a);
    if (denom == 0) return 0;

    const std::uint64_t v_min = (n + denom - 1) / denom;

    const long double root = std::sqrt(static_cast<long double>(2) * a * a - 1.0L);
    std::int64_t v_max = static_cast<std::int64_t>(std::floor(root - a));
    if (v_max <= 0) return 0;

    while ((static_cast<__int128>(v_max + 1) * (v_max + 1) + static_cast<__int128>(2) * a * (v_max + 1)) <=
           static_cast<__int128>(n)) {
        ++v_max;
    }
    while ((static_cast<__int128>(v_max) * v_max + static_cast<__int128>(2) * a * v_max) >
           static_cast<__int128>(n)) {
        --v_max;
    }

    if (v_max < static_cast<std::int64_t>(v_min)) return 0;

    factorize_int(a - 1, spf, left_factors);
    factorize_int(a + 1, spf, right_factors);
    merge_factors(left_factors, right_factors, merged_factors);

    divisors.clear();
    divisors.push_back(1);
    for (const Factor& f : merged_factors) {
        const std::size_t base_size = divisors.size();
        std::uint64_t mul = 1;
        for (int e = 1; e <= f.e; ++e) {
            mul *= static_cast<std::uint64_t>(f.p);
            for (std::size_t i = 0; i < base_size; ++i) {
                const std::uint64_t candidate = divisors[i] * mul;
                if (candidate <= static_cast<std::uint64_t>(v_max)) {
                    divisors.push_back(candidate);
                }
            }
        }
    }

    std::uint64_t count = 0;
    for (std::uint64_t v : divisors) {
        if (v < v_min) continue;

        const std::uint64_t u = n / v;
        if (((u - v) & 1ULL) != 0ULL) continue;
        if (u < v + static_cast<std::uint64_t>(2 * a)) continue;
        if (static_cast<std::uint64_t>(a) + u > static_cast<std::uint64_t>(perimeter_limit)) continue;

        ++count;
    }

    return count;
}

std::uint64_t count_barely_acute(int perimeter_limit, int threads) {
    if (perimeter_limit < 3) return 0;

    const int a_max = perimeter_limit / 3;
    const std::uint64_t a1_count = static_cast<std::uint64_t>((perimeter_limit - 1) / 2);
    const std::vector<int> spf = build_spf(a_max + 1);

    if (threads < 1) threads = 1;
    if (threads > a_max) threads = a_max;

    std::atomic<int> next_a{2};
    constexpr int kChunk = 2048;

    std::vector<std::uint64_t> partial(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            std::vector<Factor> left_factors;
            std::vector<Factor> right_factors;
            std::vector<Factor> merged_factors;
            std::vector<std::uint64_t> divisors;

            left_factors.reserve(8);
            right_factors.reserve(8);
            merged_factors.reserve(16);
            divisors.reserve(256);

            std::uint64_t local = 0;
            while (true) {
                const int start = next_a.fetch_add(kChunk, std::memory_order_relaxed);
                if (start > a_max) break;
                const int end = std::min(a_max, start + kChunk - 1);

                for (int a = start; a <= end; ++a) {
                    local += count_for_a(a,
                                         perimeter_limit,
                                         spf,
                                         left_factors,
                                         right_factors,
                                         merged_factors,
                                         divisors);
                }
            }

            partial[static_cast<std::size_t>(t)] = local;
        });
    }

    for (auto& th : pool) th.join();

    std::uint64_t total = a1_count;
    for (std::uint64_t v : partial) total += v;
    return total;
}

std::uint64_t brute_count(int perimeter_limit) {
    std::uint64_t count = 0;
    for (int a = 1; a <= perimeter_limit / 3; ++a) {
        for (int b = a; b <= (perimeter_limit - a) / 2; ++b) {
            const long long c2 = 1LL * a * a + 1LL * b * b - 1LL;
            if (c2 <= 0) continue;

            long long c = static_cast<long long>(std::sqrt(static_cast<long double>(c2)));
            while ((c + 1) * (c + 1) <= c2) ++c;
            while (c * c > c2) --c;

            if (c < b) continue;
            if (c * c != c2) continue;
            if (a + b + c > perimeter_limit) continue;
            if (a + b <= c) continue;

            ++count;
        }
    }
    return count;
}

bool validate() {
    const std::vector<int> brute_limits = {80, 120, 200, 300, 500};
    for (int limit : brute_limits) {
        const std::uint64_t slow = brute_count(limit);
        const std::uint64_t fast = count_barely_acute(limit, 1);
        if (slow != fast) {
            std::cerr << "Validation failed at perimeter=" << limit
                      << ": brute=" << slow << ", fast=" << fast << "\n";
            return false;
        }
    }

    const int thread_check_limit = 50'000;
    const std::uint64_t single = count_barely_acute(thread_check_limit, 1);

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 2;
    const int multi_threads = static_cast<int>(std::min<unsigned>(hw, 8));
    const std::uint64_t multi = count_barely_acute(thread_check_limit, multi_threads);
    if (single != multi) {
        std::cerr << "Thread consistency failed at perimeter=" << thread_check_limit
                  << ": single=" << single << ", multi=" << multi << "\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!validate()) {
        return 1;
    }

    int perimeter_limit = kDefaultLimit;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;

    if (argc > 1) perimeter_limit = std::max(3, std::atoi(argv[1]));
    if (argc > 2) threads = std::max(1, std::atoi(argv[2]));

    const std::uint64_t answer = count_barely_acute(perimeter_limit, threads);
    std::cout << answer << '\n';
    return 0;
}
