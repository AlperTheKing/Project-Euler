#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <pthread.h>
#include <unistd.h>
#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using bigint = boost::multiprecision::cpp_int;

bool extinct_for_k(int c, int d, u64 k) {
    if (k == 0) return true;

    const int s = c + d;
    std::array<u64, 321> window{};
    window[static_cast<std::size_t>(s - 1)] = k;

    int zero_count = s - 1;
    int head = 0;
    int tap = s - d;
    u64 steps = 0;
    constexpr u64 CHECK_MASK = 63ULL;

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

bool extinct_for_k_c1(int n, u64 k) {
    if (k == 0) return true;

    const int s = n + 1;
    std::array<u64, 240> window{};
    window[static_cast<std::size_t>(s - 1)] = k;

    int zero_count = s - 1;
    for (;;) {
        for (int i = 0; i < s; ++i) {
            const int j = (i + 1 == s) ? 0 : (i + 1);
            const u64 old = window[static_cast<std::size_t>(i)];
            const u64 next = (window[static_cast<std::size_t>(i)] + window[static_cast<std::size_t>(j)]) >> 1U;

            window[static_cast<std::size_t>(i)] = next;
            zero_count += static_cast<int>(old != 0 && next == 0) - static_cast<int>(old == 0 && next != 0);
        }

        if (zero_count == s) return true;
        if (zero_count == 0) return false;
    }
}

u64 compute_k(int c, int d) {
    u64 lo = 0;
    u64 hi = 1;

    while (extinct_for_k(c, d, hi)) {
        lo = hi;
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

    return lo;
}

u64 estimate_k1(int n) {
    static constexpr std::array<std::array<long double, 4>, 8> coef{{
        {{0.17656501777829126L, 0.35981185819773187L, -9.093343254260652L, 128.0591554042931L}},
        {{0.17791093179771922L, -0.02526258663290543L, 17.778074509864577L, -228.639325004742L}},
        {{0.17688473006194433L, 0.30195121097992167L, -8.649989398001244L, 159.32792700127638L}},
        {{0.1774911030303856L, 0.09144288721929997L, 9.76665768770627L, -161.46727161611716L}},
        {{0.17665123701842447L, 0.37639947984877686L, -14.382073443985455L, 264.68150347640244L}},
        {{0.17652781023196232L, 0.37645445228054497L, -10.389054848768398L, 121.34965708864206L}},
        {{0.17725060801487988L, 0.15363937288855242L, 7.49278183636027L, -167.17442628344003L}},
        {{0.17658285799047665L, 0.3729423126656912L, -11.25711787747357L, 179.7243004469314L}},
    }};

    const auto& c = coef[static_cast<std::size_t>(n & 7)];
    const long double x = static_cast<long double>(n);
    const long double est = ((c[0] * x + c[1]) * x + c[2]) * x + c[3];
    return est > 0.0L ? static_cast<u64>(std::llround(est)) : 0ULL;
}

u64 compute_k1_from_estimate(int n) {
    const u64 guess = estimate_k1(n);

    u64 lo = (guess > 2048ULL) ? (guess - 2048ULL) : 0ULL;
    u64 hi = guess + 2048ULL;

    while (lo > 0 && !extinct_for_k_c1(n, lo)) {
        hi = lo;
        lo >>= 1U;
    }
    while (extinct_for_k_c1(n, hi)) {
        lo = hi;
        hi <<= 1U;
    }

    while (lo + 1 < hi) {
        const u64 mid = lo + ((hi - lo) >> 1U);
        if (extinct_for_k_c1(n, mid)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return lo;
}

struct K1WorkerCtx {
    std::atomic<int>* next_n = nullptr;
    std::array<u64, 240>* k1 = nullptr;
};

void* k1_worker_main(void* raw_ctx) {
    auto* ctx = static_cast<K1WorkerCtx*>(raw_ctx);
    for (;;) {
        const int n = ctx->next_n->fetch_add(1, std::memory_order_relaxed);
        if (n > 239) break;
        (*ctx->k1)[static_cast<std::size_t>(n)] = compute_k1_from_estimate(n);
    }
    return nullptr;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::array<u64, 240> k1{};
    std::array<u64, 161> k_d1{};

    std::atomic<int> next_n(1);
    K1WorkerCtx worker_ctx{&next_n, &k1};

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    if (cpu_count < 1) cpu_count = 1;
    const int thread_count = static_cast<int>(std::min<long>(64, cpu_count * 2));

    std::vector<pthread_t> threads(static_cast<std::size_t>(thread_count));
    for (int t = 0; t < thread_count; ++t) {
        const int rc = pthread_create(&threads[static_cast<std::size_t>(t)], nullptr, k1_worker_main, &worker_ctx);
        assert(rc == 0);
    }
    for (int t = 0; t < thread_count; ++t) {
        const int rc = pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
    }

    for (int c = 1; c <= 10; ++c) {
        k_d1[static_cast<std::size_t>(c)] = compute_k(c, 1);
    }
    for (int c = 11; c <= 160; ++c) {
        k_d1[static_cast<std::size_t>(c)] = k1[static_cast<std::size_t>((c + 1) / 2)];
    }

    auto k_for_reduced = [&](int rc, int rd) -> u64 {
        if (rd == 1) return k_d1[static_cast<std::size_t>(rc)];
        return k1[static_cast<std::size_t>(rd + (rc - 1) / 2)];
    };

    auto G = [&](int c, int d) -> u64 {
        const int g = std::gcd(c, d);
        const int rc = c / g;
        const int rd = d / g;
        return 2 * k_for_reduced(rc, rd) + 1;
    };

    assert(G(2, 1) == 7);
    assert(G(1, 2) == 7);
    assert(G(3, 1) == 11);
    assert(G(2, 2) == 3);
    assert(G(1, 3) == 15);
    assert(k_d1[11] == compute_k(11, 1));
    assert(k_d1[57] == compute_k(57, 1));
    assert(k_d1[160] == compute_k(160, 1));

    bigint sum = 0;
    for (int c = 1; c <= 160; ++c) {
        for (int d = 1; d <= 160; ++d) {
            sum += G(c, d);
        }
    }

    std::cout << sum << '\n';
    return 0;
}
