#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <unistd.h>
#include <vector>

using u64 = std::uint64_t;

static inline bool is_pandigital_in_base(u64 x, int base) {
    const std::uint32_t all = (base == 32) ? 0xFFFFFFFFu : ((1u << base) - 1u);
    std::uint32_t mask = 0;
    while (x) {
        mask |= 1u << (x % static_cast<u64>(base));
        if (mask == all) return true;
        x /= static_cast<u64>(base);
    }
    return mask == all;
}

static inline bool is_super_pandigital(u64 x, int n) {
    for (int b = n - 1; b >= 2; --b) {
        if (!is_pandigital_in_base(x, b)) return false;
    }
    return true;
}

static std::vector<u64> k_smallest_super_pandigital(int n, int k) {
    std::vector<int> digits(n);
    std::iota(digits.begin(), digits.end(), 0);
    std::swap(digits[0], digits[1]);

    std::vector<u64> out;
    out.reserve(k);

    do {
        u64 v = 0;
        for (int d : digits) v = v * static_cast<u64>(n) + static_cast<u64>(d);

        if (is_super_pandigital(v, n)) {
            out.push_back(v);
            if (static_cast<int>(out.size()) == k) break;
        }
    } while (std::next_permutation(digits.begin(), digits.end()));

    return out;
}

static constexpr std::array<u64, 13> FACT = {
    1ULL, 1ULL, 2ULL, 6ULL, 24ULL, 120ULL, 720ULL, 5040ULL, 40320ULL, 362880ULL,
    3628800ULL, 39916800ULL, 479001600ULL};

static inline bool p8(u64 n) {
    std::uint32_t mask = 0;
    while (n > 7ULL) {
        mask |= 1u << static_cast<unsigned>(n & 7ULL);
        n >>= 3U;
        if (mask == 0xFFu) return true;
    }
    mask |= 1u << static_cast<unsigned>(n);
    return mask == 0xFFu;
}

static inline bool p11(u64 n) {
    std::uint32_t mask = 0;
    while (n > 10ULL) {
        const u64 q = n / 11ULL;
        mask |= 1u << static_cast<unsigned>(n - q * 11ULL);
        if (mask == 0x7FFu) return true;
        n = q;
    }
    mask |= 1u << static_cast<unsigned>(n);
    return mask == 0x7FFu;
}

static inline bool p_base(u64 n, u64 base, std::uint32_t full_mask) {
    std::uint32_t mask = 0;
    while (n >= base) {
        const u64 q = n / base;
        mask |= 1u << static_cast<unsigned>(n - q * base);
        if (mask == full_mask) return true;
        n = q;
    }
    mask |= 1u << static_cast<unsigned>(n);
    return mask == full_mask;
}

static inline bool is_super12_fast(u64 n) {
    if (!p8(n) || !p11(n)) return false;
    if (!p_base(n, 10ULL, 1023u)) return false;
    if (!p_base(n, 9ULL, 511u)) return false;
    if (!p_base(n, 7ULL, 127u)) return false;
    if (!p_base(n, 6ULL, 63u)) return false;
    if (!p_base(n, 5ULL, 31u)) return false;
    if (!p_base(n, 4ULL, 15u)) return false;
    if (!p_base(n, 3ULL, 7u)) return false;
    return p_base(n, 2ULL, 3u);
}

static inline void kth_perm12(u64 k, std::array<std::uint8_t, 12>& perm) {
    std::array<std::uint8_t, 12> avail{};
    for (int i = 0; i < 12; ++i) avail[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(i);
    int rem = 12;
    for (int pos = 0; pos < 12; ++pos) {
        const u64 f = FACT[static_cast<std::size_t>(11 - pos)];
        const int idx = static_cast<int>(k / f);
        k %= f;
        perm[static_cast<std::size_t>(pos)] = avail[static_cast<std::size_t>(idx)];
        for (int j = idx + 1; j < rem; ++j) {
            avail[static_cast<std::size_t>(j - 1)] = avail[static_cast<std::size_t>(j)];
        }
        --rem;
    }
}

static inline u64 to_base12_value(const std::array<std::uint8_t, 12>& perm) {
    u64 v = 0;
    for (int i = 0; i < 12; ++i) {
        v = v * 12ULL + static_cast<u64>(perm[static_cast<std::size_t>(i)]);
    }
    return v;
}

static inline void run_block(u64 start, u64 count, std::vector<u64>& out) {
    std::array<std::uint8_t, 12> perm{};
    kth_perm12(start, perm);
    for (u64 i = 0; i < count; ++i) {
        if (perm[0] != 0U) {
            const u64 v = to_base12_value(perm);
            if (is_super12_fast(v)) out.push_back(v);
        }
        if (i + 1ULL < count) {
            std::next_permutation(perm.begin(), perm.end());
        }
    }
}

struct WorkerArgs {
    u64 start = 0;
    u64 count = 0;
    std::vector<u64>* out = nullptr;
};

static void* worker_main(void* raw) {
    auto* args = static_cast<WorkerArgs*>(raw);
    run_block(args->start, args->count, *args->out);
    return nullptr;
}

static std::vector<u64> k_smallest_super_pandigital_12_parallel(int k) {
    const u64 total = FACT[12];
    const u64 chunk = FACT[9];
    u64 next = FACT[11];

    long hw = sysconf(_SC_NPROCESSORS_ONLN);
    if (hw < 1) hw = 1;
    unsigned thread_count = static_cast<unsigned>(hw);
    if (thread_count > 16U) thread_count = 16U;
    if (thread_count < 1U) thread_count = 1U;

    std::vector<u64> hits;
    hits.reserve(64);

    while (next < total) {
        const u64 remain = total - next;
        const u64 max_blocks = (remain + chunk - 1ULL) / chunk;
        const unsigned blocks =
            static_cast<unsigned>(std::min<u64>(static_cast<u64>(thread_count), max_blocks));

        std::vector<pthread_t> tids(static_cast<std::size_t>(blocks));
        std::vector<WorkerArgs> args(static_cast<std::size_t>(blocks));
        std::vector<std::vector<u64>> local(static_cast<std::size_t>(blocks));
        bool create_failed = false;
        unsigned started = 0;

        for (unsigned t = 0; t < blocks; ++t) {
            const u64 start = next + static_cast<u64>(t) * chunk;
            const u64 cnt = std::min<u64>(chunk, total - start);
            args[static_cast<std::size_t>(t)].start = start;
            args[static_cast<std::size_t>(t)].count = cnt;
            args[static_cast<std::size_t>(t)].out = &local[static_cast<std::size_t>(t)];
            if (pthread_create(&tids[static_cast<std::size_t>(t)], nullptr, worker_main,
                               &args[static_cast<std::size_t>(t)]) != 0) {
                create_failed = true;
                break;
            }
            ++started;
        }
        for (unsigned t = 0; t < started; ++t) {
            pthread_join(tids[static_cast<std::size_t>(t)], nullptr);
        }

        if (create_failed) {
            for (unsigned t = 0; t < blocks; ++t) {
                local[static_cast<std::size_t>(t)].clear();
                run_block(args[static_cast<std::size_t>(t)].start,
                          args[static_cast<std::size_t>(t)].count,
                          local[static_cast<std::size_t>(t)]);
            }
        }

        for (unsigned t = 0; t < blocks; ++t) {
            auto& vec = local[static_cast<std::size_t>(t)];
            hits.insert(hits.end(), vec.begin(), vec.end());
        }

        next += static_cast<u64>(blocks) * chunk;
        if (static_cast<int>(hits.size()) >= k) {
            std::sort(hits.begin(), hits.end());
            if (static_cast<int>(hits.size()) >= k) {
                hits.resize(static_cast<std::size_t>(k));
                return hits;
            }
        }
    }

    std::sort(hits.begin(), hits.end());
    if (static_cast<int>(hits.size()) > k) {
        hits.resize(static_cast<std::size_t>(k));
    }
    return hits;
}

int main() {
    {
        const auto v5 = k_smallest_super_pandigital(5, 1);
        if (v5.size() != 1 || v5[0] != 978ULL) {
            std::cerr << "Validation failed for smallest 5-super-pandigital\n";
            return 1;
        }
        const auto v10 = k_smallest_super_pandigital(10, 10);
        if (v10.size() != 10 || v10[0] != 1093265784ULL) {
            std::cerr << "Validation failed for smallest 10-super-pandigital\n";
            return 1;
        }
        const u64 sum10 = std::accumulate(v10.begin(), v10.end(), 0ULL);
        if (sum10 != 20319792309ULL) {
            std::cerr << "Validation failed for sum of 10 smallest 10-super-pandigital\n";
            return 1;
        }
    }

    const auto v12 = k_smallest_super_pandigital_12_parallel(10);
    const u64 ans = std::accumulate(v12.begin(), v12.end(), 0ULL);
    std::cout << ans << "\n";
    return 0;
}
