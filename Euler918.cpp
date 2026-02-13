#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;

std::pair<i64, i64> pair_value(i64 n) {
    if (n == 1) {
        return {1, 2};
    }
    if ((n & 1LL) == 0) {
        const auto [x, y] = pair_value(n >> 1);
        return {2 * x, x - 3 * y};
    }
    const auto [x, y] = pair_value(n >> 1);
    return {x - 3 * y, 2 * y};
}

i64 a_value(i64 n) {
    return pair_value(n).first;
}

i64 prefix_sum(i64 n) {
    if ((n & 1LL) == 0) {
        return 4 - a_value(n >> 1);
    }
    return 4 - 3 * a_value((n + 1) >> 1);
}

std::vector<i64> brute_sequence(int n) {
    std::vector<i64> a(n + 2, 0);
    a[1] = 1;
    for (int k = 1; 2 * k <= n; ++k) {
        a[2 * k] = 2 * a[k];
        if (2 * k + 1 <= n) {
            a[2 * k + 1] = a[k] - 3 * a[k + 1];
        }
    }
    return a;
}

void validate() {
    const std::array<i64, 10> first_ten = {1, 2, -5, 4, 17, -10, -17, 8, -47, 34};
    const std::vector<i64> a = brute_sequence(600);
    for (int i = 1; i <= 10; ++i) {
        assert(a[i] == first_ten[i - 1]);
    }

    std::vector<i64> pref(601, 0);
    for (int i = 1; i <= 600; ++i) {
        pref[i] = pref[i - 1] + a[i];
        assert(prefix_sum(i) == pref[i]);
    }

    assert(prefix_sum(10) == -13);
}

}  // namespace

int main() {
    validate();
    std::cout << prefix_sum(1'000'000'000'000LL) << '\n';
    return 0;
}
