#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while (r * r > x) --r;
    while ((r + 1) > r && (r + 1) * (r + 1) <= x) ++r;
    return r;
}

static void gen_smooth_rec(const std::vector<u64>& primes, std::size_t idx, u64 cur, u64 limit, std::vector<u64>& out) {
    if (idx == primes.size()) {
        out.push_back(cur);
        return;
    }
    const u64 p = primes[idx];
    u64 v = cur;
    while (true) {
        gen_smooth_rec(primes, idx + 1, v, limit, out);
        if (v > limit / p) break;
        v *= p;
    }
}

static std::vector<u64> gen_smooth(u64 limit) {
    // Prime factors allowed by multiplying dimensions by any k<=25: primes <= 23.
    const std::vector<u64> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    std::vector<u64> nums;
    nums.reserve(200000);
    gen_smooth_rec(primes, 0, 1, limit, nums);
    std::sort(nums.begin(), nums.end());
    nums.erase(std::unique(nums.begin(), nums.end()), nums.end());
    return nums;
}

struct MResult {
    std::vector<u64> M;  // M[k] = minimal area with exactly k variants (k up to max_k), or INF if missing.
    u64 max_filled = 0;
};

static MResult compute_M(u64 side_limit, int max_k) {
    const auto smooth = gen_smooth(side_limit);

    // Count, for each area A, how many factor pairs (x<=y) satisfy y*10 <= 11*x and x*y=A.
    // For 23-smooth A, every divisor is 23-smooth, so "manufacturable" matches "has such factor pair".
    std::unordered_map<u64, std::uint16_t> cnt;
    cnt.reserve(smooth.size() * 4);

    for (std::size_t i = 0; i < smooth.size(); ++i) {
        const u64 x = smooth[i];
        const u64 ymax = (11 * x) / 10;  // floor(1.1*x)
        auto it_end = std::upper_bound(smooth.begin() + static_cast<std::ptrdiff_t>(i), smooth.end(), ymax);
        const std::size_t j_end = static_cast<std::size_t>(it_end - smooth.begin());
        for (std::size_t j = i; j < j_end; ++j) {
            const u64 y = smooth[j];
            const u128 area128 = static_cast<u128>(x) * static_cast<u128>(y);
            if (area128 > std::numeric_limits<u64>::max()) continue;
            const u64 area = static_cast<u64>(area128);
            auto& v = cnt[area];
            if (v <= static_cast<std::uint16_t>(max_k)) ++v;  // clamp above max_k+1
        }
    }

    const u64 INF = std::numeric_limits<u64>::max();
    std::vector<u64> M(static_cast<std::size_t>(max_k + 1), INF);
    for (const auto& [area, k] : cnt) {
        if (k < 2 || k > max_k) continue;
        if (area < M[k]) M[k] = area;
    }

    u64 max_filled = 0;
    for (int k = 2; k <= max_k; ++k) {
        if (M[k] != INF && M[k] > max_filled) max_filled = M[k];
    }
    return {std::move(M), max_filled};
}

static u128 solve_sum_M(int max_k) {
    u64 side_limit = 1'000'000;  // start moderately; we will raise this if needed.

    while (true) {
        const auto res = compute_M(side_limit, max_k);
        const u64 INF = std::numeric_limits<u64>::max();

        bool all_found = true;
        for (int k = 2; k <= max_k; ++k) {
            if (res.M[k] == INF) {
                all_found = false;
                break;
            }
        }
        if (!all_found) {
            side_limit *= 10;
            continue;
        }

        // If all M(k) are found, ensure the side limit is large enough that every area <= maxM
        // has all its near-square factor pairs represented.
        const u64 maxM = res.max_filled;
        u64 s = isqrt_u64(maxM);
        if (s * s < maxM) ++s;  // ceil(sqrt(maxM))
        const u64 needed = (11 * s + 9) / 10;  // ceil(1.1*s)
        if (side_limit < needed) {
            side_limit = needed;
            continue;
        }

        // Statement check: M(3)=889200.
        assert(res.M[3] == 889200ULL);

        u128 sum = 0;
        for (int k = 2; k <= max_k; ++k) sum += static_cast<u128>(res.M[k]);
        return sum;
    }
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::cout << to_string_u128(solve_sum_M(100)) << '\n';
    return 0;
}
