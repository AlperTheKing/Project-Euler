#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
constexpr u64 kMod = 1'000'000'007ULL;
constexpr int kStates = 243;  // 3^5

int encode(const int pref, const int n0, const int n1, const int n2, const int fmod) {
    return ((((pref * 3 + n0) * 3 + n1) * 3 + n2) * 3 + fmod);
}

int F_of_state(const int id) {
    return id % 3;
}

u64 solve(const int d) {
    std::array<std::array<int, 3>, kStates> nxt{};

    for (int id = 0; id < kStates; ++id) {
        int x = id;
        const int fmod = x % 3;
        x /= 3;
        const int n2 = x % 3;
        x /= 3;
        const int n1 = x % 3;
        x /= 3;
        const int n0 = x % 3;
        x /= 3;
        const int pref = x % 3;

        for (int cls = 0; cls < 3; ++cls) {
            const int np = (pref + cls) % 3;
            std::array<int, 3> ncnt{n0, n1, n2};
            const int add = ncnt[static_cast<std::size_t>(np)];
            const int nf = (fmod + add) % 3;
            ncnt[static_cast<std::size_t>(np)] = (ncnt[static_cast<std::size_t>(np)] + 1) % 3;
            nxt[static_cast<std::size_t>(id)][static_cast<std::size_t>(cls)] =
                encode(np, ncnt[0], ncnt[1], ncnt[2], nf);
        }
    }

    std::array<u64, kStates> dp{};
    std::array<u64, kStates> ndp{};

    dp[static_cast<std::size_t>(encode(0, 1, 0, 0, 0))] = 1ULL;

    for (int pos = 1; pos <= d; ++pos) {
        ndp.fill(0ULL);
        const std::array<int, 3> mult = (pos == 1) ? std::array<int, 3>{3, 3, 3}
                                                    : std::array<int, 3>{4, 3, 3};

        for (int id = 0; id < kStates; ++id) {
            const u64 ways = dp[static_cast<std::size_t>(id)];
            if (ways == 0ULL) {
                continue;
            }
            for (int cls = 0; cls < 3; ++cls) {
                const int to = nxt[static_cast<std::size_t>(id)][static_cast<std::size_t>(cls)];
                const u64 add = (ways * static_cast<u64>(mult[static_cast<std::size_t>(cls)])) % kMod;
                u64& cell = ndp[static_cast<std::size_t>(to)];
                cell += add;
                if (cell >= kMod) {
                    cell -= kMod;
                }
            }
        }
        dp = ndp;
    }

    u64 ans = 0ULL;
    for (int id = 0; id < kStates; ++id) {
        if (F_of_state(id) == 0) {
            ans += dp[static_cast<std::size_t>(id)];
            if (ans >= kMod) {
                ans -= kMod;
            }
        }
    }
    return ans;
}

}  // namespace

int main() {
    assert(solve(2) == 30ULL);
    assert(solve(6) == 290'898ULL);

    std::cout << solve(100'000) << '\n';
    return 0;
}

