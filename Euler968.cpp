#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 MOD = 1'000'000'007U;
constexpr int VARS = 5;
constexpr int EDGES = 10;
constexpr int BITS = 31;
constexpr std::array<int, EDGES> POW3 = {
    1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683,
};
constexpr int STATE_COUNT = 59049;  // 3^10
constexpr std::array<std::array<int, 2>, EDGES> EDGE_ENDPOINTS = {{
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2},
    {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4},
}};
constexpr std::array<u32, VARS> BASES = {2U, 3U, 5U, 7U, 11U};

std::array<std::array<std::uint8_t, EDGES>, STATE_COUNT> CARRY_DIGITS{};
std::array<std::array<std::uint8_t, EDGES>, 32> PAIR_BITS{};
std::array<std::array<u32, 32>, BITS> BIT_FACTORS{};

u32 mod_pow(u32 base, u64 exp) {
    u64 result = 1U;
    u64 cur = base % MOD;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = (result * cur) % MOD;
        }
        cur = (cur * cur) % MOD;
        exp >>= 1ULL;
    }
    return static_cast<u32>(result);
}

void build_tables() {
    for (int s = 0; s < STATE_COUNT; ++s) {
        int t = s;
        for (int e = 0; e < EDGES; ++e) {
            CARRY_DIGITS[static_cast<std::size_t>(s)][static_cast<std::size_t>(e)] =
                static_cast<std::uint8_t>(t % 3);
            t /= 3;
        }
    }

    for (int mask = 0; mask < 32; ++mask) {
        for (int e = 0; e < EDGES; ++e) {
            const int u = EDGE_ENDPOINTS[static_cast<std::size_t>(e)][0];
            const int v = EDGE_ENDPOINTS[static_cast<std::size_t>(e)][1];
            const int bit_u = (mask >> u) & 1;
            const int bit_v = (mask >> v) & 1;
            PAIR_BITS[static_cast<std::size_t>(mask)][static_cast<std::size_t>(e)] =
                static_cast<std::uint8_t>(bit_u + bit_v);
        }
    }

    std::array<std::array<u32, BITS>, VARS> powers{};
    for (int i = 0; i < VARS; ++i) {
        powers[static_cast<std::size_t>(i)][0] = BASES[static_cast<std::size_t>(i)] % MOD;
        for (int k = 1; k < BITS; ++k) {
            const u64 v = powers[static_cast<std::size_t>(i)][static_cast<std::size_t>(k - 1)];
            powers[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                static_cast<u32>((v * v) % MOD);
        }
    }

    for (int k = 0; k < BITS; ++k) {
        for (int mask = 0; mask < 32; ++mask) {
            u64 value = 1U;
            for (int i = 0; i < VARS; ++i) {
                if ((mask >> i) & 1) {
                    value = (value * powers[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]) % MOD;
                }
            }
            BIT_FACTORS[static_cast<std::size_t>(k)][static_cast<std::size_t>(mask)] =
                static_cast<u32>(value);
        }
    }
}

u32 compute_p(const std::array<u32, EDGES>& bounds) {
    std::array<std::array<std::uint8_t, EDGES>, BITS> bound_bits{};
    for (int k = 0; k < BITS; ++k) {
        for (int e = 0; e < EDGES; ++e) {
            bound_bits[static_cast<std::size_t>(k)][static_cast<std::size_t>(e)] =
                static_cast<std::uint8_t>((bounds[static_cast<std::size_t>(e)] >> k) & 1U);
        }
    }

    std::vector<u32> cur(static_cast<std::size_t>(STATE_COUNT), 0U);
    std::vector<u32> nxt(static_cast<std::size_t>(STATE_COUNT), 0U);
    std::vector<int> active;
    std::vector<int> next_active;
    active.reserve(4096);
    next_active.reserve(4096);
    active.push_back(0);
    cur[0] = 1U;

    for (int k = 0; k < BITS; ++k) {
        next_active.clear();
        const auto& bits = bound_bits[static_cast<std::size_t>(k)];

        for (int state : active) {
            const u32 state_value = cur[static_cast<std::size_t>(state)];
            if (state_value == 0U) {
                continue;
            }
            cur[static_cast<std::size_t>(state)] = 0U;

            const auto& carry = CARRY_DIGITS[static_cast<std::size_t>(state)];

            for (int mask = 0; mask < 32; ++mask) {
                const auto& pair = PAIR_BITS[static_cast<std::size_t>(mask)];
                int next_state = 0;
                for (int e = 0; e < EDGES; ++e) {
                    const int t = static_cast<int>(pair[static_cast<std::size_t>(e)]) +
                                  static_cast<int>(carry[static_cast<std::size_t>(e)]) -
                                  static_cast<int>(bits[static_cast<std::size_t>(e)]);
                    const int next_carry = (t <= 0) ? 0 : ((t + 1) >> 1);
                    next_state += next_carry * POW3[static_cast<std::size_t>(e)];
                }

                const u64 add = (static_cast<u64>(state_value) *
                                 BIT_FACTORS[static_cast<std::size_t>(k)][static_cast<std::size_t>(mask)]) %
                                MOD;
                u32& dst = nxt[static_cast<std::size_t>(next_state)];
                if (dst == 0U) {
                    next_active.push_back(next_state);
                }
                const u64 sum = static_cast<u64>(dst) + add;
                dst = static_cast<u32>(sum >= MOD ? sum - MOD : sum);
            }
        }

        active.swap(next_active);
        cur.swap(nxt);
    }

    return cur[0];
}

u32 brute_p(const std::array<u32, EDGES>& bounds) {
    const int ub_a = static_cast<int>(std::min({bounds[0], bounds[1], bounds[2], bounds[3]}));
    const int ub_b = static_cast<int>(std::min({bounds[0], bounds[4], bounds[5], bounds[6]}));
    const int ub_c = static_cast<int>(std::min({bounds[1], bounds[4], bounds[7], bounds[8]}));
    const int ub_d = static_cast<int>(std::min({bounds[2], bounds[5], bounds[7], bounds[9]}));
    const int ub_e = static_cast<int>(std::min({bounds[3], bounds[6], bounds[8], bounds[9]}));

    std::vector<u32> p2(static_cast<std::size_t>(ub_a + 1), 1U);
    std::vector<u32> p3(static_cast<std::size_t>(ub_b + 1), 1U);
    std::vector<u32> p5(static_cast<std::size_t>(ub_c + 1), 1U);
    std::vector<u32> p7(static_cast<std::size_t>(ub_d + 1), 1U);
    std::vector<u32> p11(static_cast<std::size_t>(ub_e + 1), 1U);
    for (int i = 1; i <= ub_a; ++i) {
        p2[static_cast<std::size_t>(i)] = static_cast<u32>((2ULL * p2[static_cast<std::size_t>(i - 1)]) % MOD);
    }
    for (int i = 1; i <= ub_b; ++i) {
        p3[static_cast<std::size_t>(i)] = static_cast<u32>((3ULL * p3[static_cast<std::size_t>(i - 1)]) % MOD);
    }
    for (int i = 1; i <= ub_c; ++i) {
        p5[static_cast<std::size_t>(i)] = static_cast<u32>((5ULL * p5[static_cast<std::size_t>(i - 1)]) % MOD);
    }
    for (int i = 1; i <= ub_d; ++i) {
        p7[static_cast<std::size_t>(i)] = static_cast<u32>((7ULL * p7[static_cast<std::size_t>(i - 1)]) % MOD);
    }
    for (int i = 1; i <= ub_e; ++i) {
        p11[static_cast<std::size_t>(i)] = static_cast<u32>((11ULL * p11[static_cast<std::size_t>(i - 1)]) % MOD);
    }

    u32 total = 0U;
    for (int a = 0; a <= ub_a; ++a) {
        for (int b = 0; b <= ub_b; ++b) {
            if (static_cast<u32>(a + b) > bounds[0]) {
                continue;
            }
            for (int c = 0; c <= ub_c; ++c) {
                if (static_cast<u32>(a + c) > bounds[1] || static_cast<u32>(b + c) > bounds[4]) {
                    continue;
                }
                for (int d = 0; d <= ub_d; ++d) {
                    if (static_cast<u32>(a + d) > bounds[2] || static_cast<u32>(b + d) > bounds[5] ||
                        static_cast<u32>(c + d) > bounds[7]) {
                        continue;
                    }
                    for (int e = 0; e <= ub_e; ++e) {
                        if (static_cast<u32>(a + e) > bounds[3] || static_cast<u32>(b + e) > bounds[6] ||
                            static_cast<u32>(c + e) > bounds[8] || static_cast<u32>(d + e) > bounds[9]) {
                            continue;
                        }
                        u64 term = p2[static_cast<std::size_t>(a)];
                        term = (term * p3[static_cast<std::size_t>(b)]) % MOD;
                        term = (term * p5[static_cast<std::size_t>(c)]) % MOD;
                        term = (term * p7[static_cast<std::size_t>(d)]) % MOD;
                        term = (term * p11[static_cast<std::size_t>(e)]) % MOD;
                        total = static_cast<u32>(total + term);
                        if (total >= MOD) {
                            total -= MOD;
                        }
                    }
                }
            }
        }
    }

    return total;
}

bool run_checkpoints() {
    {
        std::array<u32, EDGES> bounds{};
        bounds.fill(2U);
        if (compute_p(bounds) != 7120U) {
            std::cerr << "Checkpoint failed: P(2,...,2) != 7120" << '\n';
            return false;
        }
    }
    {
        const std::array<u32, EDGES> bounds = {1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U, 10U};
        if (compute_p(bounds) != 799809376U) {
            std::cerr << "Checkpoint failed: P(1..10) mismatch" << '\n';
            return false;
        }
    }
    {
        std::mt19937 rng(968U);
        std::uniform_int_distribution<int> dist(0, 4);
        for (int it = 0; it < 20; ++it) {
            std::array<u32, EDGES> bounds{};
            for (int e = 0; e < EDGES; ++e) {
                bounds[static_cast<std::size_t>(e)] = static_cast<u32>(dist(rng));
            }
            const u32 fast = compute_p(bounds);
            const u32 brute = brute_p(bounds);
            if (fast != brute) {
                std::cerr << "Checkpoint failed on random small case" << '\n';
                return false;
            }
        }
    }
    return true;
}

u32 solve() {
    std::vector<u32> a(1000U, 0U);
    a[0] = 1U;
    a[1] = 7U;
    for (int i = 2; i < 1000; ++i) {
        const u64 term = (7ULL * a[static_cast<std::size_t>(i - 1)] +
                          static_cast<u64>(a[static_cast<std::size_t>(i - 2)]) *
                              a[static_cast<std::size_t>(i - 2)]) %
                         MOD;
        a[static_cast<std::size_t>(i)] = static_cast<u32>(term);
    }

    std::vector<u32> q(100U, 0U);
    std::atomic<int> next_idx{0};

    unsigned workers = std::thread::hardware_concurrency();
    if (workers == 0U) {
        workers = 4U;
    }
    workers = std::min(workers, 100U);

    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (unsigned t = 0; t < workers; ++t) {
        pool.emplace_back([&]() {
            while (true) {
                const int n = next_idx.fetch_add(1, std::memory_order_relaxed);
                if (n >= 100) {
                    break;
                }
                std::array<u32, EDGES> bounds{};
                for (int e = 0; e < EDGES; ++e) {
                    bounds[static_cast<std::size_t>(e)] =
                        a[static_cast<std::size_t>(10 * n + e)];
                }
                q[static_cast<std::size_t>(n)] = compute_p(bounds);
            }
        });
    }
    for (auto& th : pool) {
        th.join();
    }

    u32 total = 0U;
    for (u32 value : q) {
        total = static_cast<u32>(total + value);
        if (total >= MOD) {
            total -= MOD;
        }
    }
    return total;
}

}  // namespace

int main() {
    build_tables();
    if (!run_checkpoints()) {
        return 1;
    }
    std::cout << solve() << '\n';
    return 0;
}

