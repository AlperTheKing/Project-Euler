#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int kMaxK = 24;
constexpr i64 kScale = 1LL << 20;

std::string to_string_u128(u128 x) {
    if (x == 0U) {
        return "0";
    }
    std::string s;
    while (x > 0U) {
        const int d = static_cast<int>(x % 10U);
        s.push_back(static_cast<char>('0' + d));
        x /= 10U;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct ValueKey {
    i64 f_scaled = 0;
    std::uint8_t star = 0;
    std::array<std::int16_t, kMaxK + 1> up{};

    bool operator==(const ValueKey& other) const {
        return f_scaled == other.f_scaled && star == other.star && up == other.up;
    }
};

struct ValueHash {
    std::size_t operator()(const ValueKey& v) const noexcept {
        std::uint64_t h = 1469598103934665603ULL;
        const auto mix = [&](std::uint64_t x) {
            h ^= x;
            h *= 1099511628211ULL;
        };
        mix(static_cast<std::uint64_t>(v.f_scaled));
        mix(static_cast<std::uint64_t>(v.star));
        for (std::int16_t c : v.up) {
            mix(static_cast<std::uint16_t>(c));
        }
        return static_cast<std::size_t>(h);
    }
};

ValueKey add_values(const ValueKey& a, const ValueKey& b) {
    ValueKey out;
    out.f_scaled = a.f_scaled + b.f_scaled;
    out.star = static_cast<std::uint8_t>(a.star ^ b.star);
    for (int i = 0; i <= kMaxK; ++i) {
        out.up[static_cast<std::size_t>(i)] =
            static_cast<std::int16_t>(a.up[static_cast<std::size_t>(i)] + b.up[static_cast<std::size_t>(i)]);
    }
    return out;
}

ValueKey value_for_zero_fraction_case(int n3) {
    ValueKey out;
    int zeros_to_right = 0;
    int x = n3;
    while (x > 0) {
        const int d = x % 3;
        if (d == 0) {
            ++zeros_to_right;
        } else {
            assert(d == 2);
            out.star ^= 1U;
            if (zeros_to_right <= kMaxK) {
                ++out.up[static_cast<std::size_t>(zeros_to_right)];
            }
        }
        x /= 3;
    }
    return out;
}

std::vector<ValueKey> build_values(int limit) {
    std::vector<ValueKey> values(static_cast<std::size_t>(limit + 1));
    values[0] = ValueKey{};

    for (int n = 1; n <= limit; ++n) {
        const int q = n / 3;
        const int r = n % 3;
        ValueKey cur = values[static_cast<std::size_t>(q)];

        if (r == 1) {
            cur.f_scaled -= kScale;
        } else if (r == 2) {
            cur.star ^= 1U;
        } else {
            if (cur.f_scaled <= -2 * kScale) {
                cur.f_scaled += kScale;
            } else if (cur.f_scaled < 0) {
                cur.f_scaled /= 2;
            } else {
                assert(cur.f_scaled == 0);
                cur = value_for_zero_fraction_case(n);
            }
        }

        values[static_cast<std::size_t>(n)] = cur;
    }
    return values;
}

u128 solve_fast(int limit) {
    const std::vector<ValueKey> values = build_values(limit);

    std::unordered_map<ValueKey, u64, ValueHash> freq;
    freq.reserve(static_cast<std::size_t>(limit));
    for (int n = 1; n <= limit; ++n) {
        ++freq[values[static_cast<std::size_t>(n)]];
    }

    std::vector<ValueKey> uniq;
    uniq.reserve(freq.size());
    for (const auto& [key, _] : freq) {
        uniq.push_back(key);
    }

    std::vector<u64> f;
    f.reserve(uniq.size());
    for (const auto& key : uniq) {
        f.push_back(freq[key]);
    }

    std::unordered_map<ValueKey, u64, ValueHash> pair_count;
    pair_count.reserve(uniq.size() * 8ULL + 1ULL);

    for (std::size_t i = 0; i < uniq.size(); ++i) {
        const u64 fi = f[i];
        {
            const ValueKey s = add_values(uniq[i], uniq[i]);
            pair_count[s] += fi * (fi + 1ULL) / 2ULL;
        }
        for (std::size_t j = i + 1; j < uniq.size(); ++j) {
            const ValueKey s = add_values(uniq[i], uniq[j]);
            pair_count[s] += fi * f[j];
        }
    }

    u128 ans = 0U;
    for (const auto& [_, c] : pair_count) {
        ans += static_cast<u128>(c) * static_cast<u128>(c);
    }
    return ans;
}

}  // namespace

int main() {
    assert(solve_fast(5) == 21U);

    const u128 answer = solve_fast(100'000);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
