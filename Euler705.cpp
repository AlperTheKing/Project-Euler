#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mul_mod(const u64 a, const u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) %
                            static_cast<u128>(kMod));
}

struct DigitInfo {
    int count = 0;
    std::array<int, 4> divisors{};
    std::array<std::uint8_t, 10> has{};
};

std::array<DigitInfo, 10> build_digit_info() {
    std::array<DigitInfo, 10> info;
    for (int d = 1; d <= 9; ++d) {
        DigitInfo di;
        for (int x = 1; x <= d; ++x) {
            if (d % x == 0) {
                di.divisors[static_cast<std::size_t>(di.count)] = x;
                di.has[static_cast<std::size_t>(x)] = 1U;
                ++di.count;
            }
        }
        info[static_cast<std::size_t>(d)] = di;
    }
    return info;
}

u64 solve(const int N) {
    const auto info = build_digit_info();

    u64 total_sequences = 1ULL;
    u64 total_inversions = 0ULL;
    std::array<u64, 10> value_occurrence_sum{};

    auto process_digit = [&](const int digit) {
        if (digit == 0) {
            return;
        }

        const DigitInfo& di = info[static_cast<std::size_t>(digit)];
        const int s = di.count;

        std::array<u64, 10> gt{};
        u64 running = 0ULL;
        for (int t = 9; t >= 1; --t) {
            gt[static_cast<std::size_t>(t)] = running;
            running += value_occurrence_sum[static_cast<std::size_t>(t)];
            if (running >= kMod) {
                running -= kMod;
            }
        }

        u64 cross = 0ULL;
        for (int i = 0; i < s; ++i) {
            const int v = di.divisors[static_cast<std::size_t>(i)];
            cross += gt[static_cast<std::size_t>(v)];
            if (cross >= kMod) {
                cross -= kMod;
            }
        }

        total_inversions = (mul_mod(total_inversions, static_cast<u64>(s)) + cross) % kMod;

        std::array<u64, 10> next_occurrence{};
        for (int t = 1; t <= 9; ++t) {
            u64 v = mul_mod(value_occurrence_sum[static_cast<std::size_t>(t)], static_cast<u64>(s));
            if (di.has[static_cast<std::size_t>(t)] != 0U) {
                v += total_sequences;
                if (v >= kMod) {
                    v -= kMod;
                }
            }
            next_occurrence[static_cast<std::size_t>(t)] = v;
        }
        value_occurrence_sum = next_occurrence;
        total_sequences = mul_mod(total_sequences, static_cast<u64>(s));
    };

    if (N > 2) {
        std::vector<std::uint8_t> composite(static_cast<std::size_t>(N), 0U);
        const int root = static_cast<int>(std::sqrt(static_cast<long double>(N))) + 1;

        for (int p = 3; p <= root; p += 2) {
            if (composite[static_cast<std::size_t>(p)] != 0U) {
                continue;
            }
            if (static_cast<u64>(p) * static_cast<u64>(p) >= static_cast<u64>(N)) {
                break;
            }
            const int step = 2 * p;
            for (int x = p * p; x < N; x += step) {
                composite[static_cast<std::size_t>(x)] = 1U;
            }
        }

        process_digit(2);

        int buf[10];
        for (int p = 3; p < N; p += 2) {
            if (composite[static_cast<std::size_t>(p)] != 0U) {
                continue;
            }
            int x = p;
            int len = 0;
            while (x > 0) {
                buf[len++] = x % 10;
                x /= 10;
            }
            for (int i = len - 1; i >= 0; --i) {
                process_digit(buf[i]);
            }
        }
    }

    return total_inversions;
}

}  // namespace

int main() {
    assert(solve(20) == 3312ULL);
    assert(solve(50) == 338079744ULL);

    std::cout << solve(100'000'000) << '\n';
    return 0;
}

