#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 kRMod = 32'745'673LL;
constexpr i64 kRSeed = 6'563'116LL;
constexpr i64 kThreshold = 10'000'000LL;

i64 solve_case(i64 L, int N, int j) {
    i64 r = kRSeed;
    i64 cumulative_gap = 0;
    const i64 boundary = L - 10 - 20LL * static_cast<i64>(j - 1);

    std::vector<i64> times;
    times.reserve(static_cast<std::size_t>(N));

    for (int idx = 1; idx <= N; ++idx) {
        const i64 gap = (r % 1000) + 1;
        cumulative_gap += gap;

        const bool eastward = (r <= kThreshold);
        const i64 t = eastward ? (boundary - cumulative_gap) : (boundary + cumulative_gap);
        times.push_back(t);

        r = static_cast<i64>((static_cast<__int128>(r) * r) % kRMod);
    }

    const int kth = N - j;
    std::nth_element(times.begin(), times.begin() + kth, times.end());
    return times[kth];
}

}  // namespace

int main() {
    assert(solve_case(5000, 3, 2) == 5519);
    assert(solve_case(10'000, 11, 6) == 11'780);
    assert(solve_case(100'000, 101, 51) == 114'101);

    std::cout << solve_case(1'000'000'000, 1'000'001, 500'001) << "\n";
    return 0;
}
