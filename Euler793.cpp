#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 MOD = 50'515'093ULL;

static std::vector<u32> build_values(int n) {
    std::vector<u32> a(n);
    u64 s = 290'797ULL;
    for (int i = 0; i < n; ++i) {
        a[i] = static_cast<u32>(s);
        s = (s * s) % MOD;
    }
    std::sort(a.begin(), a.end());
    return a;
}

static u64 count_pairs_leq(const std::vector<u32>& a, u64 x) {
    const int n = static_cast<int>(a.size());
    int j = n - 1;
    u64 cnt = 0;

    for (int i = 0; i < n; ++i) {
        while (i < j && static_cast<u128>(a[i]) * static_cast<u128>(a[j]) > static_cast<u128>(x)) {
            --j;
        }
        if (j <= i) {
            break;
        }
        cnt += static_cast<u64>(j - i);
    }

    return cnt;
}

static u64 solve(int n) {
    auto a = build_values(n);

    const u64 m = static_cast<u64>(n) * static_cast<u64>(n - 1) / 2ULL;
    const u64 target = (m + 1ULL) / 2ULL;

    u64 lo = 0;
    u64 hi = static_cast<u64>(a[n - 1]) * static_cast<u64>(a[n - 2]);

    while (lo < hi) {
        const u64 mid = lo + (hi - lo) / 2ULL;
        if (count_pairs_leq(a, mid) >= target) {
            hi = mid;
        } else {
            lo = mid + 1ULL;
        }
    }

    return lo;
}

int main() {
    assert(solve(3) == 3'878'983'057'768ULL);
    assert(solve(103) == 492'700'616'748'525ULL);

    std::cout << solve(1'000'003) << '\n';
    return 0;
}
