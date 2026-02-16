#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <functional>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using u32 = std::uint32_t;

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) != 0 && (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

u64 icbrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while ((r + 1) != 0 && (r + 1) * (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r * r > x) {
        --r;
    }
    return r;
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const u64 digit = static_cast<u64>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

// D(n) = sum_{m<=n} tau(m) = sum_{i=1..n} floor(n/i).
u64 divisor_summatory(const u64 n) {
    u128 res = 0;
    for (u64 l = 1; l <= n;) {
        const u64 q = n / l;
        const u64 r = n / q;
        res += static_cast<u128>(q) * static_cast<u128>(r - l + 1);
        l = r + 1;
    }
    return static_cast<u64>(res);
}

struct Part1Group {
    u64 v;
    u64 sum_phi;
};

struct Part1WorkerCtx {
    const std::vector<Part1Group>* groups;
    unsigned worker_id;
    unsigned worker_count;
    u128 partial;
};

void* part1_worker_main(void* ptr) {
    auto* ctx = static_cast<Part1WorkerCtx*>(ptr);
    u128 local = 0;
    for (std::size_t i = ctx->worker_id; i < ctx->groups->size(); i += ctx->worker_count) {
        const Part1Group& g = (*ctx->groups)[i];
        local += static_cast<u128>(g.sum_phi) * static_cast<u128>(divisor_summatory(g.v));
    }
    ctx->partial = local;
    return nullptr;
}

unsigned choose_thread_count(std::size_t tasks) {
    if (tasks <= 1U) {
        return 1U;
    }
    long cpu_count = sysconf(_SC_NPROCESSORS_ONLN);
    unsigned threads = (cpu_count > 0) ? static_cast<unsigned>(cpu_count) : 1U;
    if (threads > 8U) {
        threads = 8U;
    }
    if (threads > tasks) {
        threads = static_cast<unsigned>(tasks);
    }
    return threads == 0U ? 1U : threads;
}

std::vector<u32> totient_sieve(const u64 n) {
    std::vector<u32> phi(static_cast<std::size_t>(n + 1), 0);
    for (u64 i = 0; i <= n; ++i) {
        phi[static_cast<std::size_t>(i)] = static_cast<u32>(i);
    }
    if (n >= 1) {
        phi[1] = 1;
    }
    for (u64 p = 2; p <= n; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (u64 m = p; m <= n; m += p) {
            const u32 cur = phi[static_cast<std::size_t>(m)];
            phi[static_cast<std::size_t>(m)] = static_cast<u32>(cur - cur / static_cast<u32>(p));
        }
    }
    return phi;
}

u128 compute_F(const u64 N) {
    const u64 limit = isqrt_u64(N);
    const u64 K = icbrt_u64(N);

    const auto phi = totient_sieve(limit);
    std::vector<u64> phi_prefix(static_cast<std::size_t>(limit + 1), 0ULL);
    for (u64 i = 1; i <= limit; ++i) {
        phi_prefix[static_cast<std::size_t>(i)] =
            phi_prefix[static_cast<std::size_t>(i - 1)] + static_cast<u64>(phi[static_cast<std::size_t>(i)]);
    }

    const u64 max_small = N / (K * K);  // ~= K
    std::vector<u32> tau(static_cast<std::size_t>(max_small + 1), 0);
    for (u64 d = 1; d <= max_small; ++d) {
        for (u64 m = d; m <= max_small; m += d) {
            ++tau[static_cast<std::size_t>(m)];
        }
    }
    std::vector<u64> D_small(static_cast<std::size_t>(max_small + 1), 0);
    for (u64 i = 1; i <= max_small; ++i) {
        D_small[static_cast<std::size_t>(i)] =
            D_small[static_cast<std::size_t>(i - 1)] + static_cast<u64>(tau[static_cast<std::size_t>(i)]);
    }

    // F(N) = sum_{t<=sqrt(N)} phi(t) * D(floor(N/t^2)).
    u128 ans = 0;

    // For t <= K, floor(N/t^2) is large; compute D via divisor_summatory.
    // Group equal quotients.
    std::vector<Part1Group> groups;
    groups.reserve(static_cast<std::size_t>(K));
    u64 t = 1;
    while (t <= K) {
        const u64 v = N / (t * t);
        u64 r = isqrt_u64(N / v);
        if (r > K) {
            r = K;
        }
        const u64 sum_phi =
            phi_prefix[static_cast<std::size_t>(r)] - phi_prefix[static_cast<std::size_t>(t - 1)];
        groups.push_back(Part1Group{v, sum_phi});
        t = r + 1;
    }

    const unsigned thread_count = choose_thread_count(groups.size());
    if (thread_count == 1U) {
        for (const Part1Group& g : groups) {
            ans += static_cast<u128>(g.sum_phi) * static_cast<u128>(divisor_summatory(g.v));
        }
    } else {
        std::vector<pthread_t> threads(thread_count);
        std::vector<Part1WorkerCtx> ctx(thread_count);
        unsigned created = 0U;
        bool failed = false;
        for (unsigned w = 0U; w < thread_count; ++w) {
            ctx[w] = Part1WorkerCtx{&groups, w, thread_count, 0};
            if (pthread_create(&threads[w], nullptr, part1_worker_main, &ctx[w]) != 0) {
                failed = true;
                break;
            }
            ++created;
        }
        for (unsigned w = 0U; w < created; ++w) {
            pthread_join(threads[w], nullptr);
        }
        if (failed) {
            for (const Part1Group& g : groups) {
                ans += static_cast<u128>(g.sum_phi) * static_cast<u128>(divisor_summatory(g.v));
            }
        } else {
            for (const Part1WorkerCtx& w : ctx) {
                ans += w.partial;
            }
        }
    }

    // For t > K, v <= max_small, so D can be tabled.
    for (u64 i = K + 1; i <= limit; ++i) {
        const u64 v = N / (i * i);
        ans += static_cast<u128>(phi[static_cast<std::size_t>(i)]) *
               static_cast<u128>(D_small[static_cast<std::size_t>(v)]);
    }

    return ans;
}

bool run_checkpoints() {
    if (compute_F(10) != static_cast<u128>(32)) {
        std::cerr << "Checkpoint failed: F(10)\n";
        return false;
    }
    if (compute_F(1000) != static_cast<u128>(12776)) {
        std::cerr << "Checkpoint failed: F(1000)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 N = 1'000'000'000'000'000ULL;
    const u128 ans = compute_F(N);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
