#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

constexpr u128 kStartX = 45;
constexpr u128 kStartY = 90;
constexpr int kMaxM = 120;

u64 splitmix64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31U);
}

struct StateKey {
    u128 x;
    u128 y;
    std::uint16_t r;
    std::uint16_t s;

    bool operator==(const StateKey& other) const {
        return x == other.x && y == other.y && r == other.r && s == other.s;
    }
};

struct StateKeyHash {
    std::size_t operator()(const StateKey& key) const {
        const u64 x_lo = static_cast<u64>(key.x);
        const u64 x_hi = static_cast<u64>(key.x >> 64U);
        const u64 y_lo = static_cast<u64>(key.y);
        const u64 y_hi = static_cast<u64>(key.y >> 64U);
        u64 h = splitmix64(x_lo);
        h ^= splitmix64(x_hi + 0x123456789abcdef0ULL);
        h ^= splitmix64(y_lo + 0x3141592653589793ULL);
        h ^= splitmix64(y_hi + 0x2718281828459045ULL);
        h ^= splitmix64((static_cast<u64>(key.r) << 16U) ^ static_cast<u64>(key.s));
        return static_cast<std::size_t>(h);
    }
};

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

struct SearchSolver {
    int R;
    int S;
    const std::vector<u128>& pow2;
    std::unordered_map<StateKey, std::uint8_t, StateKeyHash> memo;

    SearchSolver(const int r_target, const int s_target, const std::vector<u128>& pow2_ref)
        : R(r_target), S(s_target), pow2(pow2_ref) {
        memo.reserve(4096U);
    }

    bool can_reach_zero(const u128 x, const u128 y, const int r_rem, const int s_rem) const {
        const i128 min_diff = static_cast<i128>(x) * static_cast<i128>(pow2[static_cast<std::size_t>(s_rem)]) +
                              static_cast<i128>(r_rem) -
                              (static_cast<i128>(y) * static_cast<i128>(pow2[static_cast<std::size_t>(r_rem)]) +
                               static_cast<i128>(s_rem) * static_cast<i128>(pow2[static_cast<std::size_t>(r_rem)]));

        const i128 max_diff = static_cast<i128>(x) * static_cast<i128>(pow2[static_cast<std::size_t>(s_rem)]) +
                              static_cast<i128>(r_rem) * static_cast<i128>(pow2[static_cast<std::size_t>(s_rem)]) -
                              (static_cast<i128>(y) * static_cast<i128>(pow2[static_cast<std::size_t>(r_rem)]) +
                               static_cast<i128>(s_rem));

        return min_diff <= 0 && max_diff >= 0;
    }

    std::uint8_t count_paths(const u128 x, const u128 y, const int r_used, const int s_used) {
        if (r_used == R && s_used == S) {
            return (x == y) ? 1U : 0U;
        }
        if (x == y) {
            return 0U;
        }

        const int r_rem = R - r_used;
        const int s_rem = S - s_used;
        if (!can_reach_zero(x, y, r_rem, s_rem)) {
            return 0U;
        }

        const StateKey key{x, y,
                           static_cast<std::uint16_t>(r_used),
                           static_cast<std::uint16_t>(s_used)};
        const auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        std::uint8_t total = 0U;
        if (r_used < R) {
            total = static_cast<std::uint8_t>(
                std::min<int>(2, total + count_paths(x + 1U, y * 2U, r_used + 1, s_used)));
        }
        if (s_used < S && total < 2U) {
            total = static_cast<std::uint8_t>(
                std::min<int>(2, total + count_paths(x * 2U, y + 1U, r_used, s_used + 1)));
        }

        memo.emplace(key, total);
        return total;
    }

    std::string reconstruct_one_path() {
        std::string path;
        path.reserve(static_cast<std::size_t>(R + S));

        u128 x = kStartX;
        u128 y = kStartY;
        int r_used = 0;
        int s_used = 0;

        while (!(r_used == R && s_used == S)) {
            bool moved = false;

            if (r_used < R && count_paths(x + 1U, y * 2U, r_used + 1, s_used) > 0U) {
                path.push_back('r');
                x = x + 1U;
                y = y * 2U;
                ++r_used;
                moved = true;
            } else if (s_used < S && count_paths(x * 2U, y + 1U, r_used, s_used + 1) > 0U) {
                path.push_back('s');
                x = x * 2U;
                y = y + 1U;
                ++s_used;
                moved = true;
            }

            if (!moved) {
                return {};
            }
        }

        return path;
    }
};

u128 apply_path_and_final_value(const std::string& path) {
    u128 x = kStartX;
    u128 y = kStartY;
    for (std::size_t i = 0; i < path.size(); ++i) {
        if (path[i] == 'r') {
            x = x + 1U;
            y = y * 2U;
        } else {
            x = x * 2U;
            y = y + 1U;
        }
        if (i + 1U < path.size() && x == y) {
            return 0U;
        }
    }
    return (x == y) ? x : 0U;
}

struct SearchResult {
    int m = -1;
    int R = -1;
    int S = -1;
    std::uint8_t count = 0;
    u128 final_value = 0;
    std::string path;
};

std::optional<SearchResult> find_first_solution(const bool only_even_m,
                                                const std::vector<u128>& pow2) {
    for (int m = 0; m <= kMaxM; ++m) {
        if (only_even_m && (m % 2 != 0)) {
            continue;
        }

        SearchResult best{};
        bool found = false;

        for (int R = 0; R <= m; ++R) {
            const int S = m - R;

            SearchSolver solver(R, S, pow2);
            const std::uint8_t c = solver.count_paths(kStartX, kStartY, 0, 0);
            if (c == 0U) {
                continue;
            }

            const std::string path = solver.reconstruct_one_path();
            const u128 final_value = apply_path_and_final_value(path);
            if (final_value == 0U) {
                continue;
            }

            if (!found) {
                found = true;
                best.m = m;
                best.R = R;
                best.S = S;
                best.count = c;
                best.path = path;
                best.final_value = final_value;
            } else {
                best.count = 2U;
            }
        }

        if (found) {
            return best;
        }
    }
    return std::nullopt;
}

}  // namespace

int main() {
    std::vector<u128> pow2(static_cast<std::size_t>(kMaxM + 1), 1U);
    for (int i = 1; i <= kMaxM; ++i) {
        pow2[static_cast<std::size_t>(i)] = pow2[static_cast<std::size_t>(i - 1)] * 2U;
    }

    const auto first_any = find_first_solution(false, pow2);
    assert(first_any.has_value());
    assert(first_any->m == 9);
    assert(first_any->final_value == 1476U);

    const auto first_odd_length = find_first_solution(true, pow2);
    assert(first_odd_length.has_value());
    assert(first_odd_length->count == 1U);

    std::cout << to_string_u128(first_odd_length->final_value) << '\n';
    return 0;
}
