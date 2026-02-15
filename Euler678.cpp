#include <pthread.h>
#include <bit>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Entry {
    u64 n;
    int c;
    int f;
};

int detect_thread_count(std::size_t work_items) {
    long cores = ::sysconf(_SC_NPROCESSORS_ONLN);
    int threads = (cores > 0) ? static_cast<int>(cores) : 4;
    if (threads < 1) threads = 1;
    if (threads > 4) threads = 4;
    if (work_items > 0 && static_cast<std::size_t>(threads) > work_items) {
        threads = static_cast<int>(work_items);
    }
    return threads;
}

struct EWorkerTask {
    const std::vector<u64>* pw = nullptr;
    const std::vector<int>* amax_by_b = nullptr;
    const std::unordered_map<u64, u64>* perfect_mask = nullptr;
    u64 mask_filter = 0ULL;
    bool even_e = false;
    std::atomic<int>* next_b = nullptr;
    int max_b = 0;
    u64 partial = 0ULL;
};

void* e_worker_entry(void* raw) {
    auto* task = static_cast<EWorkerTask*>(raw);
    u64 local = 0ULL;

    while (true) {
        const int b = task->next_b->fetch_add(1, std::memory_order_relaxed);
        if (b > task->max_b) break;

        const int amax = (*task->amax_by_b)[static_cast<std::size_t>(b)];
        if (amax <= 0) continue;

        const u64 pb = (*task->pw)[static_cast<std::size_t>(b)];
        if (task->even_e && ((b & 1) != 0)) {
            for (int a = 2; a <= amax; a += 2) {
                const u64 s = pb + (*task->pw)[static_cast<std::size_t>(a)];
                const auto it = task->perfect_mask->find(s);
                if (it != task->perfect_mask->end()) {
                    local += static_cast<u64>(std::popcount(it->second & task->mask_filter));
                }
            }
        } else {
            for (int a = 1; a <= amax; ++a) {
                const u64 s = pb + (*task->pw)[static_cast<std::size_t>(a)];
                const auto it = task->perfect_mask->find(s);
                if (it != task->perfect_mask->end()) {
                    local += static_cast<u64>(std::popcount(it->second & task->mask_filter));
                }
            }
        }
    }

    task->partial = local;
    return nullptr;
}

u64 pow_limit(u64 base, int exp, u64 limit) {
    u128 v = 1;
    for (int i = 0; i < exp; ++i) {
        v *= static_cast<u128>(base);
        if (v > static_cast<u128>(limit)) {
            return limit + 1ULL;
        }
    }
    return static_cast<u64>(v);
}

std::vector<Entry> generate_entries(u64 n, int& max_f, int& max_c) {
    max_f = 0;
    max_c = 0;
    while ((1ULL << max_f) <= n && max_f < 62) {
        ++max_f;
    }

    std::vector<Entry> out;
    for (int f = 3; f <= max_f; ++f) {
        for (int c = 2;; ++c) {
            const u64 v = pow_limit(static_cast<u64>(c), f, n);
            if (v > n) {
                break;
            }
            out.push_back({v, c, f});
            if (c > max_c) {
                max_c = c;
            }
        }
    }
    return out;
}

std::vector<int> build_spf(int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));

    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const u64 v = static_cast<u64>(p) * static_cast<u64>(i);
            if (v > static_cast<u64>(n) || p > spf[static_cast<std::size_t>(i)]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }

    return spf;
}

std::vector<std::pair<int, int>> factorize_int(int x, const std::vector<int>& spf) {
    std::vector<std::pair<int, int>> fac;
    int v = x;
    while (v > 1) {
        const int p = spf[static_cast<std::size_t>(v)];
        int e = 0;
        do {
            v /= p;
            ++e;
        } while (v % p == 0);
        fac.push_back({p, e});
    }
    return fac;
}

u64 count_sum_two_squares_from_c_f(int c, int f, const std::vector<int>& spf) {
    const auto fac = factorize_int(c, spf);

    u64 prod = 1ULL;
    bool all_even_other = true;
    int e2 = 0;

    for (const auto& pe : fac) {
        const int p = pe.first;
        const int e = pe.second * f;

        if (p == 2) {
            e2 = e;
        } else if ((e & 1) != 0) {
            all_even_other = false;
        }

        if ((p & 3) == 3 && (e & 1) != 0) {
            return 0ULL;
        }
        if ((p & 3) == 1) {
            prod *= static_cast<u64>(e + 1);
        }
    }

    const u64 r2 = 4ULL * prod;
    const u64 zero_repr = (all_even_other && ((e2 & 1) == 0)) ? 4ULL : 0ULL;
    const u64 equal_repr = (all_even_other && e2 >= 1 && (((e2 - 1) & 1) == 0)) ? 4ULL : 0ULL;

    return (r2 - zero_repr - equal_repr) / 8ULL;
}

u64 isqrt_u128(u128 x) {
    long double ld = std::sqrt(static_cast<long double>(x));
    u64 r = static_cast<u64>(ld);

    while (static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > x) {
        --r;
    }
    return r;
}

std::vector<u64> divisors_from_factorization(const std::vector<std::pair<u64, int>>& fac) {
    std::vector<u64> divs{1ULL};
    for (const auto& pe : fac) {
        const u64 p = pe.first;
        const int e = pe.second;
        const std::size_t old = divs.size();
        u64 mul = 1ULL;
        for (int i = 1; i <= e; ++i) {
            mul *= p;
            for (std::size_t j = 0; j < old; ++j) {
                divs.push_back(divs[j] * mul);
            }
        }
    }
    return divs;
}

u64 count_sum_two_cubes_from_factorization(u64 n, const std::vector<std::pair<u64, int>>& fac) {
    const std::vector<u64> divs = divisors_from_factorization(fac);
    u64 count = 0ULL;

    for (u64 u : divs) {
        const u64 v = n / u;
        const u128 u2 = static_cast<u128>(u) * static_cast<u128>(u);

        if (u2 <= static_cast<u128>(v)) {
            continue;
        }

        const u128 t = u2 - static_cast<u128>(v);
        if (t % 3U != 0U) {
            continue;
        }
        const u128 p = t / 3U;

        if (u2 < 4U * p) {
            continue;
        }
        const u128 delta = u2 - 4U * p;

        const u64 s = isqrt_u128(delta);
        if (static_cast<u128>(s) * static_cast<u128>(s) != delta) {
            continue;
        }

        if (s >= u || ((u - s) & 1ULL) != 0ULL) {
            continue;
        }

        const u64 a = (u - s) / 2ULL;
        const u64 b = (u + s) / 2ULL;
        if (a == 0ULL || a >= b) {
            continue;
        }

        const u128 lhs = static_cast<u128>(a) * a * a + static_cast<u128>(b) * b * b;
        if (lhs == static_cast<u128>(n)) {
            ++count;
        }
    }

    return count;
}

u64 solve(u64 n) {
    int max_f = 0;
    int max_c = 0;
    const std::vector<Entry> entries = generate_entries(n, max_f, max_c);

    std::unordered_map<u64, u64> perfect_mask;
    perfect_mask.reserve(entries.size() * 2U);
    for (const auto& e : entries) {
        perfect_mask[e.n] |= (1ULL << e.f);
    }

    std::vector<u64> filter(static_cast<std::size_t>(max_f + 1), 0ULL);
    for (int e = 0; e <= max_f; ++e) {
        u64 m = 0ULL;
        for (int f = 3; f <= max_f; ++f) {
            if (e == 0 || (f % e) != 0) {
                m |= (1ULL << f);
            }
        }
        filter[static_cast<std::size_t>(e)] = m;
    }

    const std::vector<int> spf = build_spf(max_c);

    u64 ans = 0ULL;

    for (const auto& e : entries) {
        ans += count_sum_two_squares_from_c_f(e.c, e.f, spf);
    }

    std::unordered_map<u64, u64> cache3;
    cache3.reserve(1 << 16);

    for (const auto& e : entries) {
        if ((e.f % 3) == 0) {
            continue;
        }

        auto it = cache3.find(e.n);
        if (it != cache3.end()) {
            ans += it->second;
            continue;
        }

        const auto fac_c = factorize_int(e.c, spf);
        std::vector<std::pair<u64, int>> fac_n;
        fac_n.reserve(fac_c.size());
        for (const auto& pe : fac_c) {
            fac_n.push_back({static_cast<u64>(pe.first), pe.second * e.f});
        }

        const u64 val = count_sum_two_cubes_from_factorization(e.n, fac_n);
        cache3.emplace(e.n, val);
        ans += val;
    }

    for (int e = 4; e <= max_f; ++e) {
        std::vector<u64> pw;
        pw.push_back(0ULL);
        for (u64 b = 1;; ++b) {
            const u64 v = pow_limit(b, e, n);
            if (v > n) {
                break;
            }
            pw.push_back(v);
        }

        const int max_b = static_cast<int>(pw.size()) - 1;
        if (max_b < 2) {
            continue;
        }

        const bool even_e = ((e & 1) == 0);
        const u64 mask_filter = filter[static_cast<std::size_t>(e)];

        std::vector<int> amax_by_b(static_cast<std::size_t>(max_b + 1), 0);
        int amax_limit = max_b;
        for (int b = 2; b <= max_b; ++b) {
            const u64 pb = pw[static_cast<std::size_t>(b)];
            while (amax_limit > 0 && pw[static_cast<std::size_t>(amax_limit)] > n - pb) {
                --amax_limit;
            }
            int amax = amax_limit;
            if (amax >= b) amax = b - 1;
            amax_by_b[static_cast<std::size_t>(b)] = amax;
        }

        const int threads = detect_thread_count(static_cast<std::size_t>(max_b - 1));
        if (threads <= 1 || max_b < 3000) {
            for (int b = 2; b <= max_b; ++b) {
                const int amax = amax_by_b[static_cast<std::size_t>(b)];
                if (amax <= 0) continue;
                const u64 pb = pw[static_cast<std::size_t>(b)];
                if (even_e && ((b & 1) != 0)) {
                    for (int a = 2; a <= amax; a += 2) {
                        const u64 s = pb + pw[static_cast<std::size_t>(a)];
                        const auto it = perfect_mask.find(s);
                        if (it != perfect_mask.end()) {
                            ans += static_cast<u64>(std::popcount(it->second & mask_filter));
                        }
                    }
                } else {
                    for (int a = 1; a <= amax; ++a) {
                        const u64 s = pb + pw[static_cast<std::size_t>(a)];
                        const auto it = perfect_mask.find(s);
                        if (it != perfect_mask.end()) {
                            ans += static_cast<u64>(std::popcount(it->second & mask_filter));
                        }
                    }
                }
            }
        } else {
            std::vector<pthread_t> handles(static_cast<std::size_t>(threads));
            std::vector<EWorkerTask> tasks(static_cast<std::size_t>(threads));
            std::atomic<int> next_b{2};

            for (int t = 0; t < threads; ++t) {
                auto& task = tasks[static_cast<std::size_t>(t)];
                task.pw = &pw;
                task.amax_by_b = &amax_by_b;
                task.perfect_mask = &perfect_mask;
                task.mask_filter = mask_filter;
                task.even_e = even_e;
                task.next_b = &next_b;
                task.max_b = max_b;
                task.partial = 0ULL;
                pthread_create(&handles[static_cast<std::size_t>(t)], nullptr, e_worker_entry, &task);
            }
            for (int t = 0; t < threads; ++t) {
                pthread_join(handles[static_cast<std::size_t>(t)], nullptr);
                ans += tasks[static_cast<std::size_t>(t)].partial;
            }
        }
    }

    return ans;
}

}  // namespace

int main() {
    assert(solve(1'000ULL) == 7ULL);
    assert(solve(100'000ULL) == 53ULL);
    assert(solve(10'000'000ULL) == 287ULL);

    std::cout << solve(1'000'000'000'000'000'000ULL) << "\n";
    return 0;
}
