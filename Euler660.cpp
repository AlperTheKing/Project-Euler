#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <tuple>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct LengthTuple {
    int l1;
    int l2;
    int l3;
    u64 l1_low;
    u64 l1_high;
    u64 l2_low;
    u64 l2_high;
    u64 l3_low;
    u64 l3_high;
};

u64 ipow_u64(int base, int exp) {
    u64 value = 1;
    for (int i = 0; i < exp; ++i) value *= static_cast<u64>(base);
    return value;
}

u32 unique_digit_mask(u64 x, int base) {
    u32 mask = 0;
    while (true) {
        const int d = static_cast<int>(x % static_cast<u64>(base));
        const u32 bit = (1u << d);
        if (mask & bit) return 0;
        mask |= bit;
        x /= static_cast<u64>(base);
        if (x == 0) break;
    }
    return mask;
}

bool is_pandigital_triangle(u64 a, u64 b, u64 c, int base) {
    const u32 full = (1u << base) - 1u;
    const u32 ma = unique_digit_mask(a, base);
    if (ma == 0) return false;
    const u32 mb = unique_digit_mask(b, base);
    if (mb == 0 || (ma & mb)) return false;
    const u32 mc = unique_digit_mask(c, base);
    if (mc == 0 || ((ma | mb) & mc)) return false;
    return (ma | mb | mc) == full;
}

std::vector<LengthTuple> build_length_tuples(int base) {
    std::vector<LengthTuple> tuples;
    for (int l1 = 1; l1 <= base - 2; ++l1) {
        for (int l2 = l1; l2 <= base - l1 - 1; ++l2) {
            const int l3 = base - l1 - l2;
            if (l2 > l3) continue;
            if (l2 < l3 - 1) continue;
            tuples.push_back({
                l1,
                l2,
                l3,
                ipow_u64(base, l1 - 1),
                ipow_u64(base, l1) - 1,
                ipow_u64(base, l2 - 1),
                ipow_u64(base, l2) - 1,
                ipow_u64(base, l3 - 1),
                ipow_u64(base, l3) - 1,
            });
        }
    }
    return tuples;
}

u64 ceil_div_u64(u64 a, u64 b) {
    return (a + b - 1) / b;
}

u64 search_base_max_c(int base, int m_max) {
    const auto tuples = build_length_tuples(base);
    u64 best = 0;

    for (int m = 2; m <= m_max; ++m) {
        for (int n = 1; n < m; ++n) {
            if (std::gcd(m, n) != 1) continue;
            if ((m - n) % 3 == 0) continue;

            const u64 uu = static_cast<u64>(m) * static_cast<u64>(m) -
                           static_cast<u64>(n) * static_cast<u64>(n);
            const u64 vv = 2ULL * static_cast<u64>(m) * static_cast<u64>(n) +
                           static_cast<u64>(n) * static_cast<u64>(n);
            const u64 ww = static_cast<u64>(m) * static_cast<u64>(m) +
                           static_cast<u64>(m) * static_cast<u64>(n) +
                           static_cast<u64>(n) * static_cast<u64>(n);

            std::array<u64, 3> sides{uu, vv, ww};
            std::sort(sides.begin(), sides.end());
            const u64 s1 = sides[0];
            const u64 s2 = sides[1];
            const u64 s3 = sides[2];

            for (const LengthTuple& t : tuples) {
                u64 d_low = 1;
                d_low = std::max(d_low, ceil_div_u64(t.l1_low, s1));
                d_low = std::max(d_low, ceil_div_u64(t.l2_low, s2));
                d_low = std::max(d_low, ceil_div_u64(t.l3_low, s3));

                const u64 d_high = std::min(t.l1_high / s1, std::min(t.l2_high / s2, t.l3_high / s3));
                if (d_low > d_high) continue;

                for (u64 d = d_low; d <= d_high; ++d) {
                    const u64 a = s1 * d;
                    const u64 b = s2 * d;
                    const u64 c = s3 * d;
                    if (!is_pandigital_triangle(a, b, c, base)) continue;
                    if (c > best) best = c;
                }
            }
        }
    }

    return best;
}

u64 brute_base9_max_c() {
    std::array<int, 9> d{};
    for (int i = 0; i < 9; ++i) d[static_cast<std::size_t>(i)] = i;

    auto to_num = [&](int l, int r) -> u64 {
        u64 v = 0;
        for (int i = l; i < r; ++i) v = v * 9ULL + static_cast<u64>(d[static_cast<std::size_t>(i)]);
        return v;
    };

    u64 best = 0;
    do {
        for (int i = 1; i <= 7; ++i) {
            for (int j = i + 1; j <= 8; ++j) {
                if (d[0] == 0 || d[static_cast<std::size_t>(i)] == 0 || d[static_cast<std::size_t>(j)] == 0) continue;
                const u64 a0 = to_num(0, i);
                const u64 b0 = to_num(i, j);
                const u64 c0 = to_num(j, 9);

                std::array<u64, 3> s{a0, b0, c0};
                std::sort(s.begin(), s.end());
                const u64 a = s[0];
                const u64 b = s[1];
                const u64 c = s[2];

                if (a * a + b * b + a * b != c * c) continue;
                if (!is_pandigital_triangle(a, b, c, 9)) continue;
                if (c > best) best = c;
            }
        }
    } while (std::next_permutation(d.begin(), d.end()));
    return best;
}

u64 solve() {
    assert(217ULL * 217ULL + 248ULL * 248ULL + 217ULL * 248ULL == 403ULL * 403ULL);
    assert(is_pandigital_triangle(217, 248, 403, 9));
    assert(brute_base9_max_c() == 679ULL);

    constexpr int kMMax = 20'000;

    std::vector<u64> best(19, 0);
    std::atomic<int> next_base{9};

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = std::min<unsigned>(threads, 10);

    auto worker = [&]() {
        while (true) {
            const int base = next_base.fetch_add(1);
            if (base > 18) break;
            best[static_cast<std::size_t>(base)] = search_base_max_c(base, kMMax);
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    assert(best[9] == 679ULL);

    u64 sum = 0;
    for (int base = 9; base <= 18; ++base) sum += best[static_cast<std::size_t>(base)];
    return sum;
}

}  // namespace

int main() {
    std::cout << solve() << "\n";
    return 0;
}
