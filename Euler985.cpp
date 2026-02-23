#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <tuple>

namespace {

using i64 = std::int64_t;

struct SearchResult {
    i64 perimeter = std::numeric_limits<i64>::max();
    i64 a = 0;
    i64 b = 0;
    i64 c = 0;
};

long double clamp_unit(long double x) {
    if (x < -1.0L) {
        return -1.0L;
    }
    if (x > 1.0L) {
        return 1.0L;
    }
    return x;
}

int telescoping_depth(i64 a, i64 b, i64 c, int cap = 64) {
    const long double A = acosl(clamp_unit((static_cast<long double>(b) * b + static_cast<long double>(c) * c - static_cast<long double>(a) * a) /
                                           (2.0L * static_cast<long double>(b) * static_cast<long double>(c))));
    const long double B = acosl(clamp_unit((static_cast<long double>(a) * a + static_cast<long double>(c) * c - static_cast<long double>(b) * b) /
                                           (2.0L * static_cast<long double>(a) * static_cast<long double>(c))));

    const long double PI = acosl(-1.0L);
    long double x = A;
    long double y = B;
    long double z = PI - x - y;

    int depth = 0;
    while (depth < cap && x > 0.0L && y > 0.0L && z > 0.0L && x < PI / 2.0L && y < PI / 2.0L && z < PI / 2.0L) {
        ++depth;
        x = PI - 2.0L * x;
        y = PI - 2.0L * y;
        z = PI - 2.0L * z;
    }
    return depth;
}

std::pair<long double, long double> d_interval_for_target(int target) {
    const long double PI = acosl(-1.0L);
    long double lower = -std::numeric_limits<long double>::infinity();
    long double upper = std::numeric_limits<long double>::infinity();

    long double pow2 = 1.0L;
    for (int j = 0; j < target; ++j) {
        if ((j & 1) == 0) {
            lower = std::max(lower, -PI / (3.0L * pow2));
            upper = std::min(upper, PI / (6.0L * pow2));
        } else {
            lower = std::max(lower, -PI / (6.0L * pow2));
            upper = std::min(upper, PI / (3.0L * pow2));
        }
        pow2 *= 2.0L;
    }

    return {lower, upper};
}

SearchResult find_min_perimeter(int target) {
    const long double PI = acosl(-1.0L);
    const auto [lower, upper] = d_interval_for_target(target);
    const long double ratio_max = sinl(PI / 3.0L + upper) / sinl(PI / 3.0L + lower);

    i64 a_start = static_cast<i64>(floorl(1.0L / (ratio_max - 1.0L)));
    a_start = std::max<i64>(1, a_start - 2);
    while (static_cast<long double>(a_start) * ratio_max < static_cast<long double>(a_start + 1)) {
        ++a_start;
    }

    SearchResult best;

    for (i64 a = a_start;; ++a) {
        if (best.perimeter != std::numeric_limits<i64>::max() && 3 * a + 1 >= best.perimeter) {
            break;
        }

        const i64 c_max = static_cast<i64>(floorl(static_cast<long double>(a) * ratio_max)) + 2;
        for (i64 c = a + 1; c <= c_max; ++c) {
            const i64 b_min = std::max<i64>(a, c - a + 1);
            for (i64 b = b_min; b <= c; ++b) {
                const i64 p = a + b + c;
                if (best.perimeter != std::numeric_limits<i64>::max() && p >= best.perimeter) {
                    break;
                }
                if (telescoping_depth(a, b, c, target + 2) == target) {
                    best = {p, a, b, c};
                }
            }
        }
    }

    return best;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    assert(telescoping_depth(8, 9, 10, 8) == 2);

    const SearchResult sample = find_min_perimeter(2);
    assert(sample.perimeter == 10);
    assert(sample.a == 3 && sample.b == 3 && sample.c == 4);

    const SearchResult answer = find_min_perimeter(20);
    std::cout << answer.perimeter << '\n';
    return 0;
}
