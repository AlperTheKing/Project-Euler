#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 999'999'937ULL;

struct Pair {
    u64 x;
    u64 y;
};

u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(__builtin_sqrtl(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

u64 ceil_sqrt_u64(const u64 n) {
    const u64 r = isqrt_u64(n);
    return (r * r == n) ? r : (r + 1ULL);
}

Pair mul_pair(const Pair a, const Pair b, const u64 a_mod) {
    const u64 t1 = (a.x * b.x) % kMod;
    const u64 t2 = ((a.y * b.y) % kMod * a_mod) % kMod;
    const u64 real = (t1 + t2) % kMod;
    const u64 imag = (a.x * b.y + a.y * b.x) % kMod;
    return {real, imag};
}

u64 lucas_sum_mod(const u64 a, const u64 m, u64 n) {
    const u64 a_mod = a % kMod;
    Pair base{m % kMod, 1ULL};
    Pair result{1ULL, 0ULL};

    while (n > 0ULL) {
        if (n & 1ULL) {
            result = mul_pair(result, base, a_mod);
        }
        n >>= 1ULL;
        if (n > 0ULL) {
            base = mul_pair(base, base, a_mod);
        }
    }

    return (2ULL * result.x) % kMod;
}

u64 f_mod(const u64 a, const u64 n) {
    const u64 m = ceil_sqrt_u64(a);
    const bool is_square = (m * m == a);
    u64 value = lucas_sum_mod(a, m, n);
    if (!is_square) {
        value = (value + kMod - 1ULL) % kMod;
    }
    return value;
}

u128 f_exact_small(const u64 a, const u64 n) {
    const u64 m = ceil_sqrt_u64(a);
    const u64 d = m * m - a;

    if (d == 0ULL) {
        u128 p = 1;
        const u128 base = static_cast<u128>(2ULL * m);
        for (u64 i = 0ULL; i < n; ++i) {
            p *= base;
        }
        return p;
    }

    u128 s0 = 2;
    u128 s1 = static_cast<u128>(2ULL * m);
    if (n == 0ULL) {
        return 1;
    }
    if (n == 1ULL) {
        return s1 - 1;
    }

    for (u64 i = 2ULL; i <= n; ++i) {
        const u128 s = static_cast<u128>(2ULL * m) * s1 - static_cast<u128>(d) * s0;
        s0 = s1;
        s1 = s;
    }
    return s1 - 1;
}

u64 G_mod(const int limit) {
    auto sum_range = [](u64 l, u64 r) -> u64 {
        if (l > r) {
            return 0ULL;
        }
        u64 sum = 0ULL;
        u64 m = isqrt_u64(l);
        if (m * m < l) {
            ++m;
        }
        u64 sq = m * m;

        for (u64 a = l; a <= r; ++a) {
            while (sq < a) {
                ++m;
                sq = m * m;
            }
            const bool is_square = (sq == a);
            u64 value = lucas_sum_mod(a, m, a * a);
            if (!is_square) {
                value = (value + kMod - 1ULL) % kMod;
            }

            sum += value;
            if (sum >= kMod) {
                sum -= kMod;
            }
        }
        return sum;
    };

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    int thread_count = (cpu_count > 1) ? static_cast<int>(cpu_count) : 1;
    if (thread_count > 16) {
        thread_count = 16;
    }
    if (thread_count > limit) {
        thread_count = limit;
    }
    if (limit < 200'000) {
        thread_count = 1;
    }

    struct Task {
        u64 l = 0;
        u64 r = 0;
        u64 partial = 0;
    };
    auto worker = [](void* raw) -> void* {
        auto* t = static_cast<Task*>(raw);
        t->partial = 0ULL;
        if (t->l <= t->r) {
            u64 sum = 0ULL;
            u64 m = isqrt_u64(t->l);
            if (m * m < t->l) {
                ++m;
            }
            u64 sq = m * m;
            for (u64 a = t->l; a <= t->r; ++a) {
                while (sq < a) {
                    ++m;
                    sq = m * m;
                }
                const bool is_square = (sq == a);
                u64 value = lucas_sum_mod(a, m, a * a);
                if (!is_square) {
                    value = (value + kMod - 1ULL) % kMod;
                }
                sum += value;
                if (sum >= kMod) {
                    sum -= kMod;
                }
            }
            t->partial = sum;
        }
        return nullptr;
    };

    if (thread_count <= 1) {
        return sum_range(1ULL, static_cast<u64>(limit));
    }

    std::vector<pthread_t> tids(static_cast<std::size_t>(thread_count));
    std::vector<Task> tasks(static_cast<std::size_t>(thread_count));
    const u64 total = static_cast<u64>(limit);
    const u64 base = total / static_cast<u64>(thread_count);
    const u64 rem = total % static_cast<u64>(thread_count);
    u64 cur = 1ULL;
    for (int t = 0; t < thread_count; ++t) {
        const u64 len = base + (static_cast<u64>(t) < rem ? 1ULL : 0ULL);
        tasks[static_cast<std::size_t>(t)].l = cur;
        tasks[static_cast<std::size_t>(t)].r = (len == 0ULL ? 0ULL : cur + len - 1ULL);
        cur += len;
        const int rc = ::pthread_create(&tids[static_cast<std::size_t>(t)], nullptr, worker, &tasks[static_cast<std::size_t>(t)]);
        assert(rc == 0);
    }

    u64 sum = 0ULL;
    for (int t = 0; t < thread_count; ++t) {
        const int rc = ::pthread_join(tids[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
        sum += tasks[static_cast<std::size_t>(t)].partial;
        if (sum >= kMod) {
            sum -= kMod;
        }
    }
    return sum;
}

}  // namespace

int main() {
    assert(f_exact_small(5, 2) == 27);
    assert(f_exact_small(5, 5) == 3935);
    assert(f_mod(5, 2) == 27 % kMod);
    assert(f_mod(5, 5) == 3935 % kMod);
    assert(G_mod(1000) == 163'861'845ULL);

    std::cout << G_mod(5'000'000) << '\n';
    return 0;
}
