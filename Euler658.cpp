#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

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

u64 solve_case(int k, u64 n) {
    std::vector<u32> inv(static_cast<std::size_t>(k + 2), 0);
    inv[1] = 1;
    for (int i = 2; i <= k + 1; ++i) {
        inv[static_cast<std::size_t>(i)] =
            static_cast<u32>(kMod - mod_mul(kMod / static_cast<u64>(i),
                                            inv[static_cast<std::size_t>(kMod % static_cast<u64>(i))]));
    }

    const u64 exp = n + 1ULL;
    const u64 n1 = exp % kMod;

    auto geometric_sum = [&](int r) -> u64 {
        if (r == 0) return 1;
        if (r == 1) return n1;
        const u64 p = mod_pow(static_cast<u64>(r), exp);
        const u64 num = (p + kMod - 1) % kMod;
        return mod_mul(num, inv[static_cast<std::size_t>(r - 1)]);
    };

    u64 A_next = static_cast<u64>(k) % kMod;         // A_{k-1}
    u64 C_curr = static_cast<u64>(k + 1) % kMod;     // C(k+1, k)
    std::vector<u32> coeff(static_cast<std::size_t>(k), 0U);

    for (int r = k - 1; r >= 0; --r) {
        u64 A = 0;
        if (r == k - 1) {
            A = A_next;
        } else {
            const bool positive = (((k + 1 - r) & 1) == 0);
            const u64 signed_C = positive ? C_curr : (kMod - C_curr);
            A = (mod_mul(2, A_next) + signed_C + kMod - 1) % kMod;
            A_next = A;
        }
        coeff[static_cast<std::size_t>(r)] = static_cast<u32>(A);

        if (r > 0) {
            C_curr = mod_mul(C_curr, static_cast<u64>(r + 1));
            C_curr = mod_mul(C_curr, inv[static_cast<std::size_t>(k - r + 1)]);
        }
    }

    struct Task {
        int lo = 0;
        int hi = 0;
        u64 exp = 0;
        u64 n1 = 0;
        const std::vector<u32>* coeff = nullptr;
        const std::vector<u32>* inv = nullptr;
        u64 partial = 0;
    };

    auto worker = [](void* raw) -> void* {
        auto* t = static_cast<Task*>(raw);
        u64 sum = 0;
        for (int r = t->lo; r < t->hi; ++r) {
            u64 g = 0;
            if (r == 0) {
                g = 1;
            } else if (r == 1) {
                g = t->n1;
            } else {
                const u64 p = mod_pow(static_cast<u64>(r), t->exp);
                const u64 num = (p + kMod - 1) % kMod;
                g = mod_mul(num, (*t->inv)[static_cast<std::size_t>(r - 1)]);
            }
            sum += mod_mul((*t->coeff)[static_cast<std::size_t>(r)], g);
            if (sum >= kMod) {
                sum -= kMod;
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
    if (thread_count > k) {
        thread_count = k;
    }
    if (k < 200'000) {
        thread_count = 1;
    }

    if (thread_count <= 1) {
        Task task;
        task.lo = 0;
        task.hi = k;
        task.exp = exp;
        task.n1 = n1;
        task.coeff = &coeff;
        task.inv = &inv;
        worker(&task);
        return task.partial;
    }

    std::vector<pthread_t> tids(static_cast<std::size_t>(thread_count));
    std::vector<Task> tasks(static_cast<std::size_t>(thread_count));

    const int base = k / thread_count;
    const int rem = k % thread_count;
    int cur = 0;
    for (int t = 0; t < thread_count; ++t) {
        const int len = base + (t < rem ? 1 : 0);
        tasks[static_cast<std::size_t>(t)].lo = cur;
        tasks[static_cast<std::size_t>(t)].hi = cur + len;
        tasks[static_cast<std::size_t>(t)].exp = exp;
        tasks[static_cast<std::size_t>(t)].n1 = n1;
        tasks[static_cast<std::size_t>(t)].coeff = &coeff;
        tasks[static_cast<std::size_t>(t)].inv = &inv;
        cur += len;
        const int rc = ::pthread_create(&tids[static_cast<std::size_t>(t)], nullptr, worker, &tasks[static_cast<std::size_t>(t)]);
        assert(rc == 0);
    }

    u64 answer = 0;
    for (int t = 0; t < thread_count; ++t) {
        const int rc = ::pthread_join(tids[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
        answer += tasks[static_cast<std::size_t>(t)].partial;
        if (answer >= kMod) {
            answer -= kMod;
        }
    }
    return answer;
}

}  // namespace

int main() {
    assert(solve_case(4, 4) == 406);
    assert(solve_case(8, 8) == 27'902'680);
    assert(solve_case(10, 100) == 983'602'076);

    std::cout << solve_case(10'000'000, 1'000'000'000'000ULL) << "\n";
    return 0;
}
