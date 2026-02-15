#include <pthread.h>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unistd.h>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

u64 mod_pow(u64 base, int exp, int mod) {
    if (mod == 1) {
        return 0;
    }
    u64 result = 1 % mod;
    base %= static_cast<u64>(mod);
    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % static_cast<u64>(mod);
        }
        base = (base * base) % static_cast<u64>(mod);
        exp >>= 1;
    }
    return result;
}

u64 next_term(u64 x, int p, int mod) {
    if (p == 2) {
        return (x * x + 1ULL) % static_cast<u64>(mod);
    }
    return (mod_pow(x, p, mod) + 1ULL) % static_cast<u64>(mod);
}

bool hits_zero(int mod, int p) {
    if (mod == 1) {
        return true;
    }

    const u64 mod_u = static_cast<u64>(mod);
    const u64 start = 1ULL % mod_u;
    if (start == 0ULL) {
        return true;
    }

    u64 tortoise = start;
    u64 hare = next_term(start, p, mod);
    if (hare == 0ULL) {
        return true;
    }

    u64 power = 1ULL;
    u64 lam = 1ULL;
    while (tortoise != hare) {
        if (power == lam) {
            tortoise = hare;
            power <<= 1U;
            lam = 0ULL;
        }
        hare = next_term(hare, p, mod);
        if (hare == 0ULL) {
            return true;
        }
        ++lam;
    }

    return false;
}

struct Sieve {
    std::vector<int> primes;
    std::vector<int> spf;
};

Sieve build_sieve(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(p) * i;
            if (v > n || p > spf[i]) {
                break;
            }
            spf[static_cast<int>(v)] = p;
        }
    }

    return {std::move(primes), std::move(spf)};
}

int factor_list(int x, const std::vector<int>& spf, int* factors, const int limit) {
    int count = 0;
    while (x > 1) {
        const int p = spf[x];
        if (count < limit) {
            factors[count] = p;
        }
        ++count;

        while (x % p == 0) {
            x /= p;
        }
    }
    return count;
}

bool is_good_prime(int q, const std::vector<int>& spf) {
    if (q == 2) {
        return true;
    }

    int factor_buffer[32];
    const int divisor_count = factor_list(q - 1, spf, factor_buffer,
                                          static_cast<int>(std::size(factor_buffer)));
    std::sort(factor_buffer, factor_buffer + divisor_count, std::greater<int>());

    for (int i = 0; i < divisor_count; ++i) {
        const int p = factor_buffer[i];
        if (!hits_zero(q, p)) {
            return false;
        }
    }
    return true;
}

int detect_thread_count(std::size_t work_items) {
    long cores = ::sysconf(_SC_NPROCESSORS_ONLN);
    int threads = (cores > 0) ? static_cast<int>(cores) : 4;
    if (threads < 1) threads = 1;
    if (work_items == 0) return 1;
    if (static_cast<std::size_t>(threads) > work_items) {
        threads = static_cast<int>(work_items);
    }
    return threads;
}

struct PrimeWorkerTask {
    const std::vector<int>* primes = nullptr;
    const std::vector<int>* spf = nullptr;
    std::atomic<std::size_t>* next_idx = nullptr;
    std::vector<int> local_good;
};

void* prime_worker_entry(void* raw) {
    auto* task = static_cast<PrimeWorkerTask*>(raw);
    while (true) {
        const std::size_t idx = task->next_idx->fetch_add(1, std::memory_order_relaxed);
        if (idx >= task->primes->size()) break;
        const int q = (*task->primes)[idx];
        if (is_good_prime(q, *task->spf)) {
            task->local_good.push_back(q);
        }
    }
    return nullptr;
}

std::vector<int> good_primes_upto(int limit) {
    const Sieve sieve = build_sieve(limit);

    const int threads = detect_thread_count(sieve.primes.size());
    if (threads <= 1 || sieve.primes.size() < 5000U) {
        std::vector<int> good;
        good.reserve(sieve.primes.size() / 8);
        for (int q : sieve.primes) {
            if (is_good_prime(q, sieve.spf)) {
                good.push_back(q);
            }
        }
        return good;
    }

    std::atomic<std::size_t> next_idx{0};
    std::vector<pthread_t> handles(static_cast<std::size_t>(threads));
    std::vector<PrimeWorkerTask> tasks(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        auto& task = tasks[static_cast<std::size_t>(t)];
        task.primes = &sieve.primes;
        task.spf = &sieve.spf;
        task.next_idx = &next_idx;
        task.local_good.clear();
        task.local_good.reserve(sieve.primes.size() / (8U * static_cast<std::size_t>(threads)) + 8U);
        pthread_create(&handles[static_cast<std::size_t>(t)], nullptr, prime_worker_entry, &task);
    }

    std::vector<int> good;
    good.reserve(sieve.primes.size() / 8);
    for (int t = 0; t < threads; ++t) {
        pthread_join(handles[static_cast<std::size_t>(t)], nullptr);
        auto& vec = tasks[static_cast<std::size_t>(t)].local_good;
        good.insert(good.end(), vec.begin(), vec.end());
    }
    std::sort(good.begin(), good.end());
    return good;
}

u64 sum_squarefree_products_leq(const std::vector<int>& primes, i64 limit) {
    u64 total = 0;

    auto dfs = [&](auto&& self, std::size_t idx, i64 current) -> void {
        total += static_cast<u64>(current);
        for (std::size_t i = idx; i < primes.size(); ++i) {
            const i64 p = primes[i];
            if (current > limit / p) {
                break;
            }
            self(self, i + 1, current * p);
        }
    };

    dfs(dfs, 0, 1);
    return total;
}

u64 solve(i64 n) {
    const std::vector<int> good = good_primes_upto(static_cast<int>(n));

    for (int q : good) {
        if (static_cast<i64>(q) * q > n) {
            continue;
        }
        assert(!hits_zero(q * q, 2));
    }

    return sum_squarefree_products_leq(good, n);
}

void run_validations() {
    assert(solve(20) == 18);
    assert(solve(1000) == 2089);
}

}  // namespace

int main() {
    run_validations();
    constexpr i64 kN = 10'000'000;
    std::cout << solve(kN) << '\n';
    return 0;
}
