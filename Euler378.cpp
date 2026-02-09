#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr int kTargetN = 60000000;
constexpr u64 kMod = 1000000000000000000ULL;  // last 18 digits

struct Fenwick {
    std::vector<u64> bit;

    explicit Fenwick(const int n) : bit(static_cast<std::size_t>(n + 1), 0ULL) {}

    void add(int idx, const u64 value) {
        const int n = static_cast<int>(bit.size()) - 1;
        while (idx <= n) {
            bit[static_cast<std::size_t>(idx)] += value;
            idx += idx & -idx;
        }
    }

    u64 sum_prefix(int idx) const {
        u64 s = 0ULL;
        while (idx > 0) {
            s += bit[static_cast<std::size_t>(idx)];
            idx -= idx & -idx;
        }
        return s;
    }

    u64 sum_range(const int l, const int r) const {
        if (r < l) {
            return 0ULL;
        }
        return sum_prefix(r) - sum_prefix(l - 1);
    }
};

std::vector<u16> build_tau(const int limit) {
    std::vector<u32> lp(static_cast<std::size_t>(limit + 1), 0U);
    std::vector<u16> tau(static_cast<std::size_t>(limit + 1), 0U);
    std::vector<std::uint8_t> cnt(static_cast<std::size_t>(limit + 1), 0U);
    std::vector<int> primes;
    primes.reserve(4000000);

    tau[1] = 1U;
    for (int i = 2; i <= limit; ++i) {
        if (lp[static_cast<std::size_t>(i)] == 0U) {
            lp[static_cast<std::size_t>(i)] = static_cast<u32>(i);
            primes.push_back(i);
            tau[static_cast<std::size_t>(i)] = 2U;
            cnt[static_cast<std::size_t>(i)] = 1U;
        }

        for (const int p : primes) {
            const long long x = 1LL * i * p;
            if (x > limit) {
                break;
            }

            lp[static_cast<std::size_t>(x)] = static_cast<u32>(p);
            if (p == static_cast<int>(lp[static_cast<std::size_t>(i)])) {
                cnt[static_cast<std::size_t>(x)] = static_cast<std::uint8_t>(cnt[static_cast<std::size_t>(i)] + 1U);
                tau[static_cast<std::size_t>(x)] = static_cast<u16>(
                    tau[static_cast<std::size_t>(i)] / (cnt[static_cast<std::size_t>(i)] + 1U) *
                    (cnt[static_cast<std::size_t>(x)] + 1U));
                break;
            } else {
                cnt[static_cast<std::size_t>(x)] = 1U;
                tau[static_cast<std::size_t>(x)] = static_cast<u16>(tau[static_cast<std::size_t>(i)] * 2U);
            }
        }
    }

    return tau;
}

inline int d_triangle(const int n, const std::vector<u16>& tau) {
    if ((n & 1) == 0) {
        return static_cast<int>(tau[static_cast<std::size_t>(n / 2)]) * tau[static_cast<std::size_t>(n + 1)];
    }
    return static_cast<int>(tau[static_cast<std::size_t>(n)]) * tau[static_cast<std::size_t>((n + 1) / 2)];
}

u64 count_triples_mod(const int n, const std::vector<u16>& tau) {
    int max_d = 0;
    for (int i = 1; i <= n; ++i) {
        const int d = d_triangle(i, tau);
        if (d > max_d) {
            max_d = d;
        }
    }

    const int max_rank = max_d + 2;
    Fenwick count_fw(max_rank);
    Fenwick pair_fw(max_rank);

    u64 answer = 0ULL;
    for (int i = 1; i <= n; ++i) {
        const int d = d_triangle(i, tau);
        const int rank = d + 1;

        const u64 triples_ending_here = pair_fw.sum_range(rank + 1, max_rank);
        answer += triples_ending_here;
        if (answer >= kMod) {
            answer %= kMod;
        }

        const u64 greater_left = count_fw.sum_range(rank + 1, max_rank);
        pair_fw.add(rank, greater_left);
        count_fw.add(rank, 1ULL);
    }

    return answer % kMod;
}

bool run_checkpoints(const std::vector<u16>& tau) {
    if (d_triangle(7, tau) != 6) {
        std::cerr << "Checkpoint failed: dT(7)\n";
        return false;
    }
    if (count_triples_mod(20, tau) != 14ULL) {
        std::cerr << "Checkpoint failed: Tr(20)\n";
        return false;
    }
    if (count_triples_mod(100, tau) != 5772ULL) {
        std::cerr << "Checkpoint failed: Tr(100)\n";
        return false;
    }
    if (count_triples_mod(1000, tau) != 11174776ULL) {
        std::cerr << "Checkpoint failed: Tr(1000)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    const std::vector<u16> tau = build_tau(kTargetN + 1);

    if (!skip_checkpoints && !run_checkpoints(tau)) {
        return 2;
    }

    const u64 answer = count_triples_mod(kTargetN, tau);
    std::cout << answer << '\n';
    return 0;
}
