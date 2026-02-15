#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_mul(u64 a, u64 b) {
    return (a * b) % kMod;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) result = mod_mul(result, base);
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 solve_case(int alpha, u64 n) {
    std::vector<u32> inv(static_cast<std::size_t>(alpha + 1), 0);
    inv[1] = 1;
    for (int i = 2; i <= alpha; ++i) {
        inv[static_cast<std::size_t>(i)] = static_cast<u32>(
            kMod - mod_mul(kMod / static_cast<u64>(i), inv[static_cast<std::size_t>(kMod % static_cast<u64>(i))]));
    }

    u64 comb = 1;
    std::vector<u32> coeff(static_cast<std::size_t>(alpha), 0U);
    coeff[0] = 1U;
    for (int k = 0; k < alpha - 1; ++k) {
        comb = mod_mul(comb, static_cast<u64>(alpha - k));
        comb = mod_mul(comb, inv[static_cast<std::size_t>(k + 1)]);
        coeff[static_cast<std::size_t>(k + 1)] = static_cast<u32>(comb);
    }

    const u64 n1 = (n + 1ULL) % kMod;
    const u64 exp = n + 1ULL;

    struct Task {
        int lo = 0;
        int hi = 0;
        int alpha = 0;
        u64 exp = 0;
        u64 n1 = 0;
        const std::vector<u32>* coeff = nullptr;
        const std::vector<u32>* inv = nullptr;
        u64 partial = 0;
    };

    auto worker = [](void* raw) -> void* {
        auto* t = static_cast<Task*>(raw);
        u64 sum = 0;
        for (int k = t->lo; k < t->hi; ++k) {
            u64 geo = 0;
            if (k == 0) {
                geo = 1;
            } else if (k == 1) {
                geo = t->n1;
            } else {
                const u64 p = mod_pow(static_cast<u64>(k), t->exp);
                const u64 num = (p + kMod - 1) % kMod;
                geo = mod_mul(num, (*t->inv)[static_cast<std::size_t>(k - 1)]);
            }

            const u64 term = mod_mul((*t->coeff)[static_cast<std::size_t>(k)], geo);
            const bool positive = (((t->alpha - k + 1) & 1) == 0);
            if (positive) {
                sum += term;
                if (sum >= kMod) {
                    sum -= kMod;
                }
            } else {
                if (sum >= term) {
                    sum -= term;
                } else {
                    sum += kMod - term;
                }
            }
        }
        t->partial = sum;
        return nullptr;
    };

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    int thread_count = (cpu_count > 1) ? static_cast<int>(cpu_count) : 1;
    if (thread_count > 16) {
        thread_count = 16;
    }
    if (thread_count > alpha) {
        thread_count = alpha;
    }
    if (alpha < 200'000) {
        thread_count = 1;
    }

    if (thread_count <= 1) {
        Task task;
        task.lo = 0;
        task.hi = alpha;
        task.alpha = alpha;
        task.exp = exp;
        task.n1 = n1;
        task.coeff = &coeff;
        task.inv = &inv;
        worker(&task);
        return task.partial;
    }

    std::vector<pthread_t> tids(static_cast<std::size_t>(thread_count));
    std::vector<Task> tasks(static_cast<std::size_t>(thread_count));

    const int base = alpha / thread_count;
    const int rem = alpha % thread_count;
    int cur = 0;
    for (int t = 0; t < thread_count; ++t) {
        const int len = base + (t < rem ? 1 : 0);
        tasks[static_cast<std::size_t>(t)].lo = cur;
        tasks[static_cast<std::size_t>(t)].hi = cur + len;
        tasks[static_cast<std::size_t>(t)].alpha = alpha;
        tasks[static_cast<std::size_t>(t)].exp = exp;
        tasks[static_cast<std::size_t>(t)].n1 = n1;
        tasks[static_cast<std::size_t>(t)].coeff = &coeff;
        tasks[static_cast<std::size_t>(t)].inv = &inv;
        cur += len;
        const int rc =
            ::pthread_create(&tids[static_cast<std::size_t>(t)], nullptr, worker, &tasks[static_cast<std::size_t>(t)]);
        assert(rc == 0);
    }

    u64 ans = 0;
    for (int t = 0; t < thread_count; ++t) {
        const int rc = ::pthread_join(tids[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
        ans += tasks[static_cast<std::size_t>(t)].partial;
        if (ans >= kMod) {
            ans -= kMod;
        }
    }
    return ans;
}

}  // namespace

int main() {
    assert(solve_case(3, 0) == 1);
    assert(solve_case(3, 2) == 13);
    assert(solve_case(3, 4) == 79);

    std::cout << solve_case(10'000'000, 1'000'000'000'000ULL) << "\n";
    return 0;
}
