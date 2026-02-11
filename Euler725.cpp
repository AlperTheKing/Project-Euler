#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;

constexpr u64 kMod = 10'000'000'000'000'000ULL;

u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= kMod) {
        a %= kMod;
    }
    return a;
}

u64 S(const int n) {
    u64 answer = 0;

    for (int s = 0; s <= 9; ++s) {
        const int target = 2 * s;

        std::array<std::array<u64, 2>, 19> cnt{};
        std::array<std::array<u64, 2>, 19> sum{};

        for (int d = 1; d <= 9; ++d) {
            if (d > target) {
                break;
            }
            const int has = (d == s) ? 1 : 0;
            cnt[static_cast<std::size_t>(d)][static_cast<std::size_t>(has)] =
                add_mod(cnt[static_cast<std::size_t>(d)][static_cast<std::size_t>(has)], 1);
            sum[static_cast<std::size_t>(d)][static_cast<std::size_t>(has)] =
                add_mod(sum[static_cast<std::size_t>(d)][static_cast<std::size_t>(has)], static_cast<u64>(d));
        }

        for (int len = 1; len <= n; ++len) {
            answer = add_mod(answer, sum[static_cast<std::size_t>(target)][1]);
            if (len == n) {
                break;
            }

            std::array<std::array<u64, 2>, 19> next_cnt{};
            std::array<std::array<u64, 2>, 19> next_sum{};

            for (int current_sum = 0; current_sum <= target; ++current_sum) {
                for (int has = 0; has <= 1; ++has) {
                    const u64 c = cnt[static_cast<std::size_t>(current_sum)][static_cast<std::size_t>(has)];
                    if (c == 0) {
                        continue;
                    }
                    const u64 sv = sum[static_cast<std::size_t>(current_sum)][static_cast<std::size_t>(has)];

                    for (int d = 0; d <= 9; ++d) {
                        const int ns = current_sum + d;
                        if (ns > target) {
                            break;
                        }
                        const int nh = (has != 0 || d == s) ? 1 : 0;

                        next_cnt[static_cast<std::size_t>(ns)][static_cast<std::size_t>(nh)] =
                            add_mod(next_cnt[static_cast<std::size_t>(ns)][static_cast<std::size_t>(nh)], c);

                        const u64 appended = (sv * 10ULL + c * static_cast<u64>(d)) % kMod;
                        next_sum[static_cast<std::size_t>(ns)][static_cast<std::size_t>(nh)] =
                            add_mod(next_sum[static_cast<std::size_t>(ns)][static_cast<std::size_t>(nh)], appended);
                    }
                }
            }

            cnt = next_cnt;
            sum = next_sum;
        }
    }

    return answer;
}

}  // namespace

int main() {
    assert(S(3) == 63'270ULL);
    assert(S(7) == 85'499'991'450ULL);

    std::cout << S(2020) << '\n';
    return 0;
}
