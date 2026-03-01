#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include <pthread.h>
#include <unistd.h>
#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using bigint = boost::multiprecision::cpp_int;

struct ReducedPair {
    int c = 0;
    int d = 0;
    u64 g = 0;
};

bool extinct_for_k(int c, int d, u64 k) {
    if (k == 0) return true;

    const int s = c + d;
    std::array<u64, 321> window{};
    window[static_cast<std::size_t>(s - 1)] = k;

    int zero_count = s - 1;
    int head = 0;
    int tap = s - d;
    u64 steps = 0;
    constexpr u64 CHECK_MASK = 255ULL;

    for (;;) {
        const u64 old = window[static_cast<std::size_t>(head)];
        const u64 next = (window[static_cast<std::size_t>(head)] + window[static_cast<std::size_t>(tap)]) >> 1U;

        window[static_cast<std::size_t>(head)] = next;
        zero_count += static_cast<int>(old != 0 && next == 0) - static_cast<int>(old == 0 && next != 0);
        if (++head == s) head = 0;
        if (++tap == s) tap = 0;

        ++steps;
        if ((steps & CHECK_MASK) == 0) {
            if (zero_count == s) return true;
            if (zero_count == 0) return false;
        }
    }
}

u64 compute_g_coprime(int c, int d) {
    u64 lo = 0;
    u64 hi = 1;

    while (extinct_for_k(c, d, hi)) {
        lo = hi;
        if (hi > (std::numeric_limits<u64>::max() >> 1U)) break;
        hi <<= 1U;
    }

    while (lo + 1 < hi) {
        const u64 mid = lo + ((hi - lo) >> 1U);
        if (extinct_for_k(c, d, mid)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 2 * lo + 1;
}

u64 G_from_lookup(int c, int d, const std::array<std::array<u64, 161>, 161>& lookup) {
    const int g = std::gcd(c, d);
    return lookup[static_cast<std::size_t>(c / g)][static_cast<std::size_t>(d / g)];
}

struct WorkerCtx {
    std::atomic<int>* next_index = nullptr;
    std::vector<ReducedPair>* pairs = nullptr;
};

void* worker_main(void* raw_ctx) {
    auto* ctx = static_cast<WorkerCtx*>(raw_ctx);
    for (;;) {
        const int idx = ctx->next_index->fetch_add(1, std::memory_order_relaxed);
        if (idx >= static_cast<int>(ctx->pairs->size())) break;
        auto& p = (*ctx->pairs)[static_cast<std::size_t>(idx)];
        p.g = compute_g_coprime(p.c, p.d);
    }
    return nullptr;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::vector<ReducedPair> pairs;
    pairs.reserve(16'000);
    for (int c = 1; c <= 160; ++c) {
        for (int d = 1; d <= 160; ++d) {
            if (std::gcd(c, d) == 1) pairs.push_back(ReducedPair{c, d, 0});
        }
    }

    std::atomic<int> next_index(0);
    WorkerCtx ctx{&next_index, &pairs};

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    if (cpu_count < 1) cpu_count = 1;
    const int thread_count = static_cast<int>(cpu_count);

    std::vector<pthread_t> threads(static_cast<std::size_t>(thread_count));
    for (int t = 0; t < thread_count; ++t) {
        const int rc = pthread_create(&threads[static_cast<std::size_t>(t)], nullptr, worker_main, &ctx);
        assert(rc == 0);
    }
    for (int t = 0; t < thread_count; ++t) {
        const int rc = pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
    }

    std::array<std::array<u64, 161>, 161> lookup{};
    for (const auto& p : pairs) {
        lookup[static_cast<std::size_t>(p.c)][static_cast<std::size_t>(p.d)] = p.g;
    }

    assert(extinct_for_k(2, 1, 3));
    assert(!extinct_for_k(2, 1, 4));

    assert(G_from_lookup(2, 1, lookup) == 7);
    assert(G_from_lookup(1, 2, lookup) == 7);
    assert(G_from_lookup(3, 1, lookup) == 11);
    assert(G_from_lookup(2, 2, lookup) == 3);
    assert(G_from_lookup(1, 3, lookup) == 15);

    bigint sum = 0;
    for (int c = 1; c <= 160; ++c) {
        for (int d = 1; d <= 160; ++d) {
            sum += G_from_lookup(c, d, lookup);
        }
    }

    std::cout << sum << '\n';
    return 0;
}
