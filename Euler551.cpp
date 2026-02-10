#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

static constexpr int DMAX = 12;                 // H fits in <= 12 decimal digits for this problem.
static constexpr int CMAX = 9 * DMAX;           // max possible digit sum of H.
static constexpr u32 P = 1'000'000U;            // 10^6, keeps carry per step <= 1.
static constexpr int DS6_MAX = 54;              // max digit sum for 6 digits.
static constexpr int RMAX = DS6_MAX + CMAX;     // max remainder after crossing P.

static int digit_sum_u64(u64 x) {
    int s = 0;
    while (x) {
        s += static_cast<int>(x % 10);
        x /= 10;
    }
    return s;
}

struct Tables {
    std::vector<u8> ds6;  // digit sum for 0..P-1
    // For fixed c in [0..CMAX] and starting x in [0..RMAX], jump to the next carry:
    // after steps[c][x] single-step updates, the low part overflows and becomes rem[c][x],
    // and the high part increases by 1.
    std::array<std::array<u16, RMAX + 1>, CMAX + 1> rem{};
    std::array<std::array<u32, RMAX + 1>, CMAX + 1> steps{};

    Tables() { build(); }

    void build() {
        ds6.assign(P, 0);
        for (u32 i = 1; i < P; ++i) {
            ds6[i] = static_cast<u8>(ds6[i / 10] + (i % 10));
        }

        std::vector<u32> st(P);
        std::vector<u16> rm(P);

        for (int c = 0; c <= CMAX; ++c) {
            for (int i = static_cast<int>(P) - 1; i >= 0; --i) {
                const u32 x = static_cast<u32>(i);
                const u32 y = x + static_cast<u32>(ds6[x]) + static_cast<u32>(c);
                if (y >= P) {
                    st[x] = 1;
                    rm[x] = static_cast<u16>(y - P);
                } else {
                    st[x] = st[y] + 1;
                    rm[x] = rm[y];
                }
            }
            for (int x = 0; x <= RMAX; ++x) {
                steps[c][x] = st[static_cast<u32>(x)];
                rem[c][x] = rm[static_cast<u32>(x)];
            }
        }
    }
};

struct Trans {
    std::array<u16, RMAX + 1> next{};
    std::array<u64, RMAX + 1> cost{};  // number of single-step updates
};

static Trans identity_trans() {
    Trans t;
    for (int r = 0; r <= RMAX; ++r) {
        t.next[r] = static_cast<u16>(r);
        t.cost[r] = 0;
    }
    return t;
}

static Trans compose(const Trans& a, const Trans& b) {
    // Apply a then b.
    Trans out;
    for (int r = 0; r <= RMAX; ++r) {
        const u16 mid = a.next[r];
        out.next[r] = b.next[mid];
        out.cost[r] = a.cost[r] + b.cost[mid];
    }
    return out;
}

struct BlockBuilder {
    const Tables& tab;
    std::array<u64, DMAX + 1> pow10{};
    std::array<std::array<Trans, CMAX + 1>, DMAX + 1> block{};

    explicit BlockBuilder(const Tables& t) : tab(t) { build(); }

    Trans single(int digit_sum_h) const {
        Trans tr;
        for (int r = 0; r <= RMAX; ++r) {
            tr.next[r] = tab.rem[digit_sum_h][r];
            tr.cost[r] = tab.steps[digit_sum_h][r];
        }
        return tr;
    }

    void build() {
        pow10[0] = 1;
        for (int d = 1; d <= DMAX; ++d) pow10[d] = pow10[d - 1] * 10ULL;

        for (int ps = 0; ps <= CMAX; ++ps) block[0][ps] = single(ps);

        for (int d = 1; d <= DMAX; ++d) {
            for (int ps = 0; ps + 9 * d <= CMAX; ++ps) {
                Trans tr = block[d - 1][ps];
                for (int dig = 1; dig <= 9; ++dig) tr = compose(tr, block[d - 1][ps + dig]);
                block[d][ps] = tr;
            }
        }
    }

    Trans range_trans(u64 l, u64 r) const {
        // Compose macro-steps for H in [l, r), in increasing order.
        Trans res = identity_trans();
        u64 cur = l;
        while (cur < r) {
            int best_d = 0;
            for (int d = DMAX; d >= 0; --d) {
                const u64 p = pow10[d];
                if (p == 0) continue;
                if (cur % p == 0 && cur + p <= r) {
                    best_d = d;
                    break;
                }
            }
            const u64 p = pow10[best_d];
            const u64 hi = cur / p;
            const int ps = digit_sum_u64(hi);
            res = compose(res, block[best_d][ps]);
            cur += p;
        }
        return res;
    }
};

static u64 solve_a(u64 n) {
    // a_0=1, a_1=1, and for n>=1: a_{n+1}=a_n+sumdigits(a_n).
    if (n == 0) return 1;
    if (n == 1) return 1;

    static const Tables tables;
    static const BlockBuilder builder(tables);

    u64 remaining_steps = n - 1;  // number of single-step updates from a_1.

    u64 H = 0;
    u32 low = 1;

    // Find maximum number of low-part carry events we can perform within remaining_steps.
    // Each carry event processes steps[digit_sum(H)][low] single steps and increases H by 1.
    auto cost_after_carries = [&](u64 carry_cnt, u16 start_low, u16& out_low) -> u64 {
        const Trans tr = builder.range_trans(H, H + carry_cnt);
        out_low = tr.next[start_low];
        return tr.cost[start_low];
    };

    u64 lo = 0, hi = 1;
    while (true) {
        u16 dummy = 0;
        const u64 cost = cost_after_carries(hi, static_cast<u16>(low), dummy);
        if (cost > remaining_steps) break;
        lo = hi;
        hi <<= 1ULL;
    }

    while (lo + 1 < hi) {
        const u64 mid = lo + (hi - lo) / 2;
        u16 dummy = 0;
        const u64 cost = cost_after_carries(mid, static_cast<u16>(low), dummy);
        if (cost <= remaining_steps) lo = mid;
        else hi = mid;
    }

    u16 low_after = 0;
    const u64 used = cost_after_carries(lo, static_cast<u16>(low), low_after);
    H += lo;
    remaining_steps -= used;
    low = static_cast<u32>(low_after);

    // Finish remaining steps without crossing the next carry boundary.
    const int c = digit_sum_u64(H);
    for (u64 i = 0; i < remaining_steps; ++i) {
        low += static_cast<u32>(tables.ds6[low]) + static_cast<u32>(c);
        // By construction we should not overflow.
        assert(low < P);
    }

    return H * static_cast<u64>(P) + static_cast<u64>(low);
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Checkpoint from the problem statement.
    assert(solve_a(1'000'000ULL) == 31054319ULL);

    std::cout << solve_a(1'000'000'000'000'000ULL) << '\n';
    return 0;
}

