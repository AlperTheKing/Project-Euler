#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 ceil_div_sqrt(u64 d, u64 D) {
    long double approx = static_cast<long double>(d) / std::sqrt(static_cast<long double>(D));
    u64 g = static_cast<u64>(approx);
    if (g == 0) {
        g = 1;
    }

    const u128 dd = static_cast<u128>(d) * static_cast<u128>(d);
    while (static_cast<u128>(g) * static_cast<u128>(g) * static_cast<u128>(D) < dd) {
        ++g;
    }
    while (g > 1 && static_cast<u128>(g - 1) * static_cast<u128>(g - 1) * static_cast<u128>(D) >= dd) {
        --g;
    }
    return g;
}

struct Entry {
    u32 value;
    u64 count;
};

struct State {
    u64 s = 0;
    u64 C = 0;
    u64 D = 0;
    u64 z = 0;
    std::vector<Entry> groups;

    State() = default;

    State(u64 s_, u64 C_, u64 D_, const std::vector<Entry>& input_groups)
        : s(s_), C(C_), D(D_) {
        groups = input_groups;
        normalize();
    }

    void normalize() {
        std::vector<Entry> cleaned;
        cleaned.reserve(groups.size());
        for (const auto& e : groups) {
            if (e.value > 0 && e.count > 0) {
                cleaned.push_back(e);
            }
        }

        std::sort(cleaned.begin(), cleaned.end(),
                  [](const Entry& a, const Entry& b) { return a.value < b.value; });

        groups.clear();
        groups.reserve(cleaned.size());
        for (const auto& e : cleaned) {
            if (!groups.empty() && groups.back().value == e.value) {
                groups.back().count += e.count;
            } else {
                groups.push_back(e);
            }
        }

        u64 pos_count = 0;
        for (const auto& e : groups) {
            pos_count += e.count;
        }
        assert(pos_count <= s);
        z = s - pos_count;
    }

    u64 prefix_smallest(u64 b) const {
        if (b <= z) {
            return 0ULL;
        }

        u64 remaining = b - z;
        u64 sum = 0;

        for (const auto& e : groups) {
            if (remaining == 0) {
                break;
            }
            const u64 take = std::min(remaining, e.count);
            sum += take * static_cast<u64>(e.value);
            if (sum > C) {
                return C + 1ULL;
            }
            remaining -= take;
        }

        if (remaining > 0) {
            return C + 1ULL;
        }
        return sum;
    }

    std::pair<State, u64> advance(u64 d, u64 b, u64 g, u64 cost) const {
        if (b == 0) {
            return {State(s + d, C, D, {{static_cast<u32>(C), 1ULL}}), C};
        }

        const u64 z_sel = std::min(b, z);
        u64 rem = b - z_sel;

        std::vector<Entry> selected;
        selected.reserve(groups.size());

        for (const auto& e : groups) {
            if (rem == 0) {
                break;
            }
            const u64 take = std::min(rem, e.count);
            if (take > 0) {
                selected.push_back({e.value, take});
                rem -= take;
            }
        }
        assert(rem == 0);

        const u64 c = C - cost;

        std::vector<Entry> next_entries;
        next_entries.reserve(selected.size() + 2);

        if (c > 0) {
            next_entries.push_back({static_cast<u32>(c), 1ULL});
        }
        if (z_sel > 0) {
            next_entries.push_back({static_cast<u32>(g), z_sel});
        }
        for (const auto& e : selected) {
            next_entries.push_back({static_cast<u32>(static_cast<u64>(e.value) + g), e.count});
        }

        return {State(s + d, C, D, next_entries), c};
    }
};

u128 compute_T(u64 N, u64 C, u64 D) {
    State st(1ULL, C, D, {{static_cast<u32>(C), 1ULL}});
    u64 current_c = C;
    u128 total = C;

    while (st.s < N) {
        const u64 s = st.s;
        const u64 lo = (s > 2ULL * C) ? (s - 2ULL * C) : 1ULL;

        u64 found_d = s;
        u64 found_b = 0;
        u64 found_g = 1;
        u64 found_cost = 0;

        u64 g = ceil_div_sqrt(lo, D);

        for (u64 d = lo; d <= s; ++d) {
            const u128 dd = static_cast<u128>(d) * static_cast<u128>(d);
            while (static_cast<u128>(g) * static_cast<u128>(g) * static_cast<u128>(D) < dd) {
                ++g;
            }

            const u64 b = (s - d + 1ULL) / 2ULL;
            if (b == 0) {
                found_d = d;
                found_b = 0;
                found_g = g;
                found_cost = 0;
                break;
            }

            const u64 ps = st.prefix_smallest(b);
            if (ps > C) {
                continue;
            }

            const u128 cost128 = static_cast<u128>(b) * static_cast<u128>(g) + static_cast<u128>(ps);
            if (cost128 <= C) {
                found_d = d;
                found_b = b;
                found_g = g;
                found_cost = static_cast<u64>(cost128);
                break;
            }
        }

        const u64 next_s = s + found_d;
        if (next_s > N) {
            const u64 cnt = N - s;
            total += static_cast<u128>(cnt) * static_cast<u128>(current_c);
            total += static_cast<u128>(cnt) * static_cast<u128>(cnt + 1ULL) / 2ULL;
            break;
        }

        const u64 cnt = next_s - s - 1ULL;
        if (cnt > 0) {
            total += static_cast<u128>(cnt) * static_cast<u128>(current_c);
            total += static_cast<u128>(cnt) * static_cast<u128>(cnt + 1ULL) / 2ULL;
        }

        auto [next_state, next_c] = st.advance(found_d, found_b, found_g, found_cost);
        total += next_c;
        st = std::move(next_state);
        current_c = next_c;
    }

    return total;
}

void run_validations() {
    assert(compute_T(30ULL, 3ULL, 3ULL) == 190ULL);
    assert(compute_T(50ULL, 3ULL, 31ULL) == 385ULL);
    assert(compute_T(1'000ULL, 101ULL, 101ULL) == 142'427ULL);
}

}  // namespace

int main() {
    run_validations();

    constexpr u64 kN = 10'000'000'000'000'000ULL;
    u128 sum = 0;

    u64 pow10 = 1;
    for (int k = 1; k <= 6; ++k) {
        pow10 *= 10ULL;
        const u64 C = pow10 + 1ULL;
        sum += compute_T(kN, C, C);
    }

    const u64 answer = static_cast<u64>(sum % 1'000'000'000ULL);
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}
