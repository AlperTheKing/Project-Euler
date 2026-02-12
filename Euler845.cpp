#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr int kMaxLen = 20;
static constexpr int kMaxSum = 9 * kMaxLen;

struct DigitPrimeDP {
    std::array<bool, kMaxSum + 1> is_prime{};
    std::array<std::array<u128, kMaxSum + 1>, kMaxLen + 1> ways{};
    std::array<std::array<u128, kMaxSum + 1>, kMaxLen + 1> good{};

    DigitPrimeDP() {
        sieve_primes();
        build_ways();
        build_good();
    }

    void sieve_primes() {
        is_prime.fill(true);
        is_prime[0] = false;
        is_prime[1] = false;
        for (int p = 2; p * p <= kMaxSum; ++p) {
            if (!is_prime[p]) continue;
            for (int x = p * p; x <= kMaxSum; x += p) is_prime[x] = false;
        }
    }

    void build_ways() {
        for (auto& row : ways) row.fill(0);
        ways[0][0] = 1;
        for (int len = 1; len <= kMaxLen; ++len) {
            for (int s = 0; s <= 9 * (len - 1); ++s) {
                u128 cur = ways[len - 1][s];
                if (cur == 0) continue;
                for (int d = 0; d <= 9; ++d) ways[len][s + d] += cur;
            }
        }
    }

    void build_good() {
        for (auto& row : good) row.fill(0);
        for (int rem = 0; rem <= kMaxLen; ++rem) {
            for (int pref = 0; pref <= kMaxSum; ++pref) {
                u128 cnt = 0;
                for (int tail = 0; tail <= 9 * rem; ++tail) {
                    if (pref + tail <= kMaxSum && is_prime[pref + tail]) cnt += ways[rem][tail];
                }
                good[rem][pref] = cnt;
            }
        }
    }

    u128 count_upto(u64 x) const {
        std::string s = std::to_string(x);
        int sum = 0;
        u128 ans = 0;

        for (int i = 0; i < static_cast<int>(s.size()); ++i) {
            int d = s[i] - '0';
            int rem = static_cast<int>(s.size()) - i - 1;
            for (int dig = 0; dig < d; ++dig) ans += good[rem][sum + dig];
            sum += d;
        }
        if (is_prime[sum]) ++ans;
        return ans;
    }
};

static u64 nth_value(u64 n, const DigitPrimeDP& dp) {
    u128 target = n;
    u64 lo = 1;
    u64 hi = 1;
    while (dp.count_upto(hi) < target) hi <<= 1;

    while (lo < hi) {
        u64 mid = lo + (hi - lo) / 2;
        if (dp.count_upto(mid) >= target) hi = mid;
        else lo = mid + 1;
    }
    return lo;
}

int main() {
    DigitPrimeDP dp;
    assert(nth_value(61, dp) == 157ULL);
    assert(nth_value(100'000'000ULL, dp) == 403'539'364ULL);

    std::cout << nth_value(10'000'000'000'000'000ULL, dp) << '\n';
    return 0;
}
