#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using u32 = uint32_t;
using u64 = uint64_t;

constexpr u32 kMod = 135707531u;
constexpr u32 kDefaultN = 12344321u;

u32 mod_pow(u32 base, u32 exp) {
    u64 res = 1;
    u64 cur = base % kMod;
    while (exp > 0) {
        if (exp & 1u) res = (res * cur) % kMod;
        cur = (cur * cur) % kMod;
        exp >>= 1u;
    }
    return static_cast<u32>(res);
}

u32 compute_F(u32 n) {
    if (n < 2) return 0;

    std::vector<u32> inv(n + 1, 0);
    inv[1] = 1;
    for (u32 i = 2; i <= n; ++i) {
        inv[i] = kMod - static_cast<u64>(kMod / i) * inv[kMod % i] % kMod;
    }

    u32 fact_n = 1;
    for (u32 i = 2; i <= n; ++i) {
        fact_n = static_cast<u64>(fact_n) * i % kMod;
    }

    std::vector<u32> invfact(n + 1, 1);
    for (u32 i = 1; i <= n; ++i) {
        invfact[i] = static_cast<u64>(invfact[i - 1]) * inv[i] % kMod;
    }

    u32 term_z = static_cast<u64>(n % kMod) * ((n - 1) % kMod) % kMod;
    term_z = static_cast<u64>(term_z) * mod_pow(n - 2, n - 1) % kMod;

    if (n < 4) return term_z;

    const u32 total_terms = n - 3;  // m in [2, n-2]
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;
    if (total_terms < 200000) threads = 1;
    threads = std::min<unsigned>(threads, 8);
    threads = std::min<unsigned>(threads, static_cast<unsigned>(total_terms));

    std::vector<u32> partial(threads, 0);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    const u32 chunk = (total_terms + static_cast<u32>(threads) - 1) / static_cast<u32>(threads);

    auto worker = [&](unsigned idx, u32 start, u32 end) {
        u32 local = 0;
        for (u32 m = start; m < end; ++m) {
            u32 term = static_cast<u64>(fact_n) * invfact[m] % kMod;
            term = static_cast<u64>(term) * inv[n - m] % kMod;
            term = static_cast<u64>(term) * mod_pow(m - 1, m) % kMod;
            local += term;
            if (local >= kMod) local -= kMod;
        }
        partial[idx] = local;
    };

    u32 start = 2;
    for (unsigned t = 0; t < threads; ++t) {
        u32 end = std::min(start + chunk, n - 1);  // end is exclusive, last m is n-2.
        pool.emplace_back(worker, t, start, end);
        start = end;
    }
    for (auto& th : pool) th.join();

    u32 sum_c0 = 0;
    for (u32 v : partial) {
        sum_c0 += v;
        if (sum_c0 >= kMod) sum_c0 -= kMod;
    }

    u32 total = term_z + sum_c0;
    if (total >= kMod) total -= kMod;
    return total;
}

void run_checks() {
    struct Check { u32 n; u32 expected; };
    const Check checks[] = {
        {3, 6},
        {8, 16276736},
        {100, 84326147},
    };
    for (const auto& chk : checks) {
        u32 got = compute_F(chk.n);
        if (got != chk.expected) {
            std::cerr << "Validation failed for n=" << chk.n << ": got " << got
                      << ", expected " << chk.expected << '\n';
            std::exit(1);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    run_checks();
    u32 n = kDefaultN;
    if (argc > 1) {
        n = static_cast<u32>(std::strtoul(argv[1], nullptr, 10));
    }
    std::cout << compute_F(n) << '\n';
    return 0;
}
