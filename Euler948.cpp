#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u128 binom_u128(u64 n, u64 k) {
    if (k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }

    u128 result = 1;
    for (u64 i = 1; i <= k; ++i) {
        result = (result * static_cast<u128>(n - k + i)) / static_cast<u128>(i);
    }
    return result;
}

u128 catalan_u128(u64 k) {
    return binom_u128(2 * k, k) / static_cast<u128>(k + 1);
}

u128 pow2_u128(u64 n) {
    return static_cast<u128>(1) << n;
}

u128 F_formula(u64 n) {
    if (n % 2ULL == 1ULL) {
        const u64 k = (n - 1ULL) / 2ULL;
        return pow2_u128(n) - 2ULL * binom_u128(2ULL * k, k);
    }

    const u64 k = n / 2ULL;
    return pow2_u128(n) - 2ULL * binom_u128(2ULL * k - 1ULL, k - 1ULL) - catalan_u128(k - 1ULL);
}

bool bit_is_r(u64 mask, int i) {
    return ((mask >> i) & 1ULL) != 0ULL;
}

u64 F_bruteforce(int n) {
    u64 count = 0;
    const u64 total = 1ULL << n;

    for (u64 mask = 0; mask < total; ++mask) {
        std::vector<std::vector<unsigned char>> wl(n, std::vector<unsigned char>(n, 0));
        std::vector<std::vector<unsigned char>> wr(n, std::vector<unsigned char>(n, 0));

        for (int i = 0; i < n; ++i) {
            const bool r = bit_is_r(mask, i);
            wl[i][i] = static_cast<unsigned char>(!r);
            wr[i][i] = static_cast<unsigned char>(r);
        }

        for (int len = 2; len <= n; ++len) {
            for (int i = 0; i + len - 1 < n; ++i) {
                const int j = i + len - 1;

                unsigned char left_win = 0;
                for (int t = i + 1; t <= j; ++t) {
                    if (!wr[t][j]) {
                        left_win = 1;
                        break;
                    }
                }
                wl[i][j] = left_win;

                unsigned char right_win = 0;
                for (int t = i; t < j; ++t) {
                    if (!wl[i][t]) {
                        right_win = 1;
                        break;
                    }
                }
                wr[i][j] = right_win;
            }
        }

        if (wl[0][n - 1] && wr[0][n - 1]) {
            ++count;
        }
    }

    return count;
}

void run_validations() {
    assert(F_formula(3) == 4);
    assert(F_formula(8) == 181);

    for (int n = 1; n <= 12; ++n) {
        assert(F_formula(static_cast<u64>(n)) == F_bruteforce(n));
    }
}

void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }

    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        std::cout << *it;
    }
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kN = 60;
    print_u128(F_formula(kN));
    std::cout << '\n';
    return 0;
}
