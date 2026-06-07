#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 MOD = 1'000'000'007ULL;

struct AndCandidate {
    u64 value = 0;
    u64 upper_bound = 0;
};

int bit_count(const int limit) {
    int bits = 0;
    while ((1 << bits) <= limit) {
        ++bits;
    }
    return bits;
}

int highest_bit(const int value) {
    int bit = 0;
    while ((1 << (bit + 1)) <= value) {
        ++bit;
    }
    return bit;
}

int parity_sign(int value) {
    int parity = 0;
    while (value > 0) {
        parity ^= 1;
        value &= value - 1;
    }
    return parity == 0 ? 1 : -1;
}

AndCandidate parity_and_candidate(const int n) {
    const int bits = bit_count(n);
    std::vector<int> count(static_cast<std::size_t>(bits), 0);
    std::vector<int> discrepancy(static_cast<std::size_t>(bits), 0);

    for (int value = 1; value <= n; ++value) {
        const int sign = parity_sign(value);
        for (int bit = 0; bit < bits; ++bit) {
            if (((value >> bit) & 1) != 0) {
                ++count[static_cast<std::size_t>(bit)];
                discrepancy[static_cast<std::size_t>(bit)] += sign;
            }
        }
    }

    AndCandidate result;
    for (int bit = 0; bit < bits; ++bit) {
        const u64 weight = 1ULL << bit;
        const u64 m = static_cast<u64>(count[static_cast<std::size_t>(bit)]);
        const int s = discrepancy[static_cast<std::size_t>(bit)];
        result.upper_bound += weight * ((m * m) / 4);
        result.value += weight * ((m * m - static_cast<u64>(s * s)) / 4);
    }

    return result;
}

u64 max_and_value(const int n) {
    const AndCandidate candidate = parity_and_candidate(n);
    assert(candidate.value == candidate.upper_bound);
    return candidate.value;
}

u64 brute_max_and(const int n) {
    u64 best = 0;
    for (u64 mask = 0; mask < (1ULL << n); ++mask) {
        u64 current = 0;
        for (int a = 1; a <= n; ++a) {
            for (int b = a + 1; b <= n; ++b) {
                const bool side_a = ((mask >> (a - 1)) & 1ULL) != 0;
                const bool side_b = ((mask >> (b - 1)) & 1ULL) != 0;
                if (side_a != side_b) {
                    current += static_cast<u64>(a & b);
                }
            }
        }
        best = std::max(best, current);
    }
    return best;
}

struct Edge {
    int weight;
    int from;
    int to;

    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

u64 max_xor_sum(const int n) {
    std::vector<int> squares(static_cast<std::size_t>(n + 1), 0);
    for (int value = 1; value <= n; ++value) {
        squares[static_cast<std::size_t>(value)] = value * value;
    }

    std::vector<Edge> edges;
    edges.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1));
    for (int from = 1; from <= n; ++from) {
        for (int to = 1; to <= n; ++to) {
            if (from != to) {
                edges.push_back({squares[static_cast<std::size_t>(from)] ^
                                     squares[static_cast<std::size_t>(to)],
                                 from, to});
            }
        }
    }

    std::sort(edges.begin(), edges.end());

    std::vector<u64> best(static_cast<std::size_t>(n + 1), 0);
    std::vector<std::pair<int, u64>> pending;
    u64 answer = 0;

    for (std::size_t index = 0; index < edges.size();) {
        const int weight = edges[index].weight;
        pending.clear();

        while (index < edges.size() && edges[index].weight == weight) {
            const Edge& edge = edges[index];
            const u64 candidate = best[static_cast<std::size_t>(edge.from)] + static_cast<u64>(weight);
            pending.push_back({edge.to, candidate});
            answer = std::max(answer, candidate);
            ++index;
        }

        for (const auto& update : pending) {
            u64& slot = best[static_cast<std::size_t>(update.first)];
            slot = std::max(slot, update.second);
        }
    }

    return answer;
}

u64 count_unreachable_nim(const int n) {
    if (n <= 1) {
        return 0;
    }

    const int top_bit = highest_bit(n - 1);
    u64 total = 0;
    for (int target_bit = 0; target_bit <= top_bit; ++target_bit) {
        std::array<u64, 8> dp{};
        dp[7] = 1;

        for (int bit = top_bit; bit >= 0; --bit) {
            std::array<u64, 8> next{};
            const int limit_bit = ((n - 1) >> bit) & 1;

            for (int mask = 0; mask < 8; ++mask) {
                if (dp[static_cast<std::size_t>(mask)] == 0) {
                    continue;
                }
                for (int a = 0; a <= 1; ++a) {
                    for (int b = 0; b <= 1; ++b) {
                        for (int c = 0; c <= 1; ++c) {
                            if (bit > target_bit && ((a ^ b ^ c) != 0)) {
                                continue;
                            }
                            if (bit == target_bit && (a == 0 || b == 0 || c == 0)) {
                                continue;
                            }

                            const std::array<int, 3> chosen{a, b, c};
                            bool valid = true;
                            int next_mask = 0;
                            for (int index = 0; index < 3; ++index) {
                                if (((mask >> index) & 1) == 0) {
                                    continue;
                                }
                                if (chosen[static_cast<std::size_t>(index)] > limit_bit) {
                                    valid = false;
                                    break;
                                }
                                if (chosen[static_cast<std::size_t>(index)] == limit_bit) {
                                    next_mask |= 1 << index;
                                }
                            }

                            if (valid) {
                                next[static_cast<std::size_t>(next_mask)] +=
                                    dp[static_cast<std::size_t>(mask)];
                            }
                        }
                    }
                }
            }

            dp = next;
        }

        for (const u64 count : dp) {
            total += count;
        }
    }

    return total;
}

u64 brute_unreachable_nim(const int n) {
    u64 total = 0;
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
            for (int c = 0; c < n; ++c) {
                const int nim_sum = a ^ b ^ c;
                if (nim_sum == 0) {
                    continue;
                }
                const int bit = highest_bit(nim_sum);
                if (((a >> bit) & 1) != 0 && ((b >> bit) & 1) != 0 && ((c >> bit) & 1) != 0) {
                    ++total;
                }
            }
        }
    }
    return total;
}

u64 mul_mod(const u64 lhs, const u64 rhs) {
    return static_cast<u64>((static_cast<u128>(lhs) * rhs) % MOD);
}

u64 meta_value(const int index, const u64 first, const u64 second, const u64 third) {
    if (index == 0) {
        return first % MOD;
    }
    if (index == 1) {
        return second % MOD;
    }
    if (index == 2) {
        return third % MOD;
    }

    u64 a = first % MOD;
    u64 b = second % MOD;
    u64 c = third % MOD;
    for (int current = 3; current <= index; ++current) {
        const u64 next = mul_mod(mul_mod(c, b), a);
        a = b;
        b = c;
        c = next;
    }

    return c;
}

void run_checkpoints() {
    assert(brute_max_and(10) == 50);
    assert(max_and_value(10) == 50);
    const AndCandidate target_and = parity_and_candidate(1000);
    assert(target_and.value == target_and.upper_bound);
    assert(max_xor_sum(4) == 71);
    assert(max_xor_sum(10) == 702);

    for (int n = 1; n <= 20; ++n) {
        assert(count_unreachable_nim(n) == brute_unreachable_nim(n));
    }
    assert(count_unreachable_nim(10) == 123);
}

}  // namespace

int main() {
    run_checkpoints();

    const u64 max_and = max_and_value(1000);
    const u64 max_xor = max_xor_sum(1000);
    const u64 unreachable_nim = count_unreachable_nim(1000);
    assert(meta_value(4, max_and, max_xor, unreachable_nim) == 457'587'170ULL);

    std::cout << meta_value(1000, max_and, max_xor, unreachable_nim) << '\n';
    return 0;
}
