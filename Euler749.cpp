#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_set>

namespace {

using u64 = std::uint64_t;

struct Solver {
    static constexpr int kMaxDigits = 16;
    static constexpr int kMaxExponent = 53;
    static constexpr u64 kLimit = 9'999'999'999'999'999ULL;

    std::array<u64, 10> pow{};
    std::array<int, 10> counts{};
    std::unordered_set<u64> found;

    bool matches_counts(u64 n, const int len) const {
        if (n == 0ULL) {
            return false;
        }
        std::array<int, 10> freq{};
        int digits = 0;
        while (n > 0ULL) {
            ++freq[static_cast<int>(n % 10ULL)];
            n /= 10ULL;
            ++digits;
        }
        return digits == len && freq == counts;
    }

    void check_candidate(const u64 n, const int len) {
        if (n == 0ULL || n > kLimit) {
            return;
        }
        if (matches_counts(n, len)) {
            found.insert(n);
        }
    }

    void dfs(const int digit, const int rem, const int len, const u64 sum) {
        if (sum > kLimit + 1ULL) {
            return;
        }
        if (digit == 10) {
            if (rem != 0) {
                return;
            }
            if (counts[0] == len) {
                return;
            }
            if (sum > 0ULL) {
                check_candidate(sum - 1ULL, len);
            }
            check_candidate(sum + 1ULL, len);
            return;
        }

        for (int c = 0; c <= rem; ++c) {
            counts[digit] = c;
            const u64 add = static_cast<u64>(c) * pow[digit];
            if (sum + add < sum) {
                continue;
            }
            dfs(digit + 1, rem - c, len, sum + add);
        }
        counts[digit] = 0;
    }

    u64 solve(const int max_digits) {
        found.clear();
        found.reserve(1U << 16U);

        for (int k = 1; k <= kMaxExponent; ++k) {
            pow[0] = 0ULL;
            pow[1] = 1ULL;
            for (int d = 2; d <= 9; ++d) {
                u64 v = 1ULL;
                for (int i = 0; i < k; ++i) {
                    if (v > (kLimit + 1ULL) / static_cast<u64>(d)) {
                        v = kLimit + 2ULL;
                        break;
                    }
                    v *= static_cast<u64>(d);
                }
                pow[d] = v;
            }

            for (int len = 1; len <= max_digits; ++len) {
                counts.fill(0);
                dfs(0, len, len, 0ULL);
            }
        }

        u64 total = 0ULL;
        for (const u64 n : found) {
            total += n;
        }
        return total;
    }
};

}  // namespace

int main() {
    Solver solver;
    assert(solver.solve(2) == 110ULL);
    assert(solver.solve(6) == 2'562'701ULL);

    std::cout << solver.solve(16) << '\n';
    return 0;
}
