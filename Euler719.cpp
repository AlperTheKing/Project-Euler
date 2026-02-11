#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

bool is_s_number(const i64 root) {
    const u64 square = static_cast<u64>(root) * static_cast<u64>(root);
    const std::string s = std::to_string(square);
    const int len = static_cast<int>(s.size());

    std::vector<int> digits(static_cast<std::size_t>(len));
    for (int i = 0; i < len; ++i) {
        digits[static_cast<std::size_t>(i)] = s[static_cast<std::size_t>(i)] - '0';
    }

    std::vector<i64> suffix_value(static_cast<std::size_t>(len + 1), 0);
    std::vector<int> suffix_digit_sum(static_cast<std::size_t>(len + 1), 0);

    u64 place = 1;
    for (int i = len - 1; i >= 0; --i) {
        suffix_value[static_cast<std::size_t>(i)] =
            static_cast<i64>(digits[static_cast<std::size_t>(i)] * place +
                             suffix_value[static_cast<std::size_t>(i + 1)]);
        place *= 10;
        suffix_digit_sum[static_cast<std::size_t>(i)] =
            suffix_digit_sum[static_cast<std::size_t>(i + 1)] + digits[static_cast<std::size_t>(i)];
    }

    const auto dfs = [&](const auto& self, const int pos, const i64 sum) -> bool {
        if (pos == len) {
            return sum == root;
        }

        i64 part = 0;
        for (int end = pos; end < len; ++end) {
            part = part * 10 + digits[static_cast<std::size_t>(end)];
            const i64 next_sum = sum + part;
            if (next_sum > root) {
                break;
            }

            const i64 min_possible = next_sum + suffix_digit_sum[static_cast<std::size_t>(end + 1)];
            if (min_possible > root) {
                continue;
            }
            const i64 max_possible = next_sum + suffix_value[static_cast<std::size_t>(end + 1)];
            if (max_possible < root) {
                continue;
            }

            if (self(self, end + 1, next_sum)) {
                return true;
            }
        }
        return false;
    };

    i64 first_part = 0;
    for (int end = 0; end <= len - 2; ++end) {
        first_part = first_part * 10 + digits[static_cast<std::size_t>(end)];
        if (first_part > root) {
            break;
        }

        const i64 min_possible = first_part + suffix_digit_sum[static_cast<std::size_t>(end + 1)];
        if (min_possible > root) {
            continue;
        }
        const i64 max_possible = first_part + suffix_value[static_cast<std::size_t>(end + 1)];
        if (max_possible < root) {
            continue;
        }

        if (dfs(dfs, end + 1, first_part)) {
            return true;
        }
    }

    return false;
}

i64 T(const i64 limit) {
    const i64 max_root = static_cast<i64>(1'000'000);
    i64 sum = 0;

    for (i64 root = 2; root <= max_root; ++root) {
        const int mod9 = static_cast<int>(root % 9);
        if (mod9 != 0 && mod9 != 1) {
            continue;
        }
        if (root * root > limit) {
            break;
        }
        if (is_s_number(root)) {
            sum += root * root;
        }
    }
    return sum;
}

}  // namespace

int main() {
    assert(T(10'000) == 41'333);

    std::cout << T(1'000'000'000'000LL) << '\n';
    return 0;
}
