#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

struct Tri {
    std::uint32_t a;
    std::uint32_t b;
    std::uint32_t c;
};

struct FlatTriSet {
    std::vector<Tri> table;
    u64 mask;
    std::size_t used = 0;

    explicit FlatTriSet(std::size_t pow2_capacity)
        : table(pow2_capacity, Tri{0, 0, 0}), mask(pow2_capacity - 1) {}

    static u64 hash_tri(std::uint32_t a, std::uint32_t b, std::uint32_t c) {
        u64 x = (static_cast<u64>(a) << 1) ^ (static_cast<u64>(b) << 23) ^ (static_cast<u64>(c) << 41);
        x ^= x >> 33;
        x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33;
        x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= x >> 33;
        return x;
    }

    bool insert(std::uint32_t a, std::uint32_t b, std::uint32_t c) {
        u64 idx = hash_tri(a, b, c) & mask;
        while (true) {
            Tri& cur = table[idx];
            if (cur.a == 0) {
                cur = Tri{a, b, c};
                ++used;
                return true;
            }
            if (cur.a == a && cur.b == b && cur.c == c) {
                return false;
            }
            idx = (idx + 1) & mask;
        }
    }
};

std::size_t choose_capacity(int perimeter_limit) {
    std::size_t target = static_cast<std::size_t>(3.1L * perimeter_limit);
    if (target < 4096) {
        target = 4096;
    }
    std::size_t cap = 1;
    while (cap < target) {
        cap <<= 1U;
    }
    return cap;
}

i64 fortunate_sum_fast(int perimeter_limit) {
    FlatTriSet seen(choose_capacity(perimeter_limit));

    const i64 lim = 120LL * perimeter_limit;
    const int max_m = static_cast<int>(std::sqrt(static_cast<long double>(lim))) + 2;

    for (int m = 1; m <= max_m; ++m) {
        const i64 m2 = 1LL * m * m;
        if (m2 >= lim) {
            break;
        }
        const int nmax = static_cast<int>(std::sqrt(static_cast<long double>((lim - m2) / 15LL)));
        for (int n = 1; n <= nmax; ++n) {
            if (std::gcd(m, n) != 1) {
                continue;
            }

            const i64 n2 = 1LL * n * n;
            const i64 x0 = m2 - 15LL * n2;
            if (x0 <= 0) {
                continue;
            }
            const i64 y0 = 2LL * m * n;
            const i64 z0 = m2 + 15LL * n2;

            const i64 g = std::gcd(x0, std::gcd(y0, z0));
            const i64 x = x0 / g;
            const i64 y = y0 / g;
            const i64 z = z0 / g;

            for (int s : {-1, 1}) {
                const i64 u_num = x + s * y;
                if (u_num <= 0) {
                    continue;
                }

                const i64 perimeter_num = x + z + (s == 1 ? 5LL : 3LL) * y;
                if (perimeter_num <= 0) {
                    continue;
                }
                const i64 tmax = (4LL * perimeter_limit) / perimeter_num;
                if (tmax <= 0) {
                    continue;
                }

                const i64 t_need_u = 4LL / std::gcd(4LL, u_num);
                const i64 t_need_z = 4LL / std::gcd(4LL, z);
                const i64 t0 = std::lcm(t_need_u, t_need_z);

                for (i64 t = t0; t <= tmax; t += t0) {
                    i64 a = (t * u_num) / 4LL;
                    i64 b = t * y;
                    i64 c = (t * z) / 4LL;
                    if (a > b) {
                        std::swap(a, b);
                    }
                    if (b > c) {
                        std::swap(b, c);
                    }
                    if (a > b) {
                        std::swap(a, b);
                    }
                    if (a + b + c > perimeter_limit) {
                        continue;
                    }
                    seen.insert(static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(b),
                                static_cast<std::uint32_t>(c));
                }
            }
        }
    }

    i64 answer = 0;
    for (const Tri& t : seen.table) {
        if (t.a == 0) {
            continue;
        }
        answer += static_cast<i64>(t.a) + static_cast<i64>(t.b) + static_cast<i64>(t.c);
    }
    return answer;
}

bool fortunate_angle(i64 opposite, i64 side1, i64 side2) {
    const i64 num = side1 * side1 + side2 * side2 - opposite * opposite;
    return 2LL * std::llabs(num) == side1 * side2;
}

bool is_fortunate_triangle(i64 a, i64 b, i64 c) {
    return fortunate_angle(a, b, c) || fortunate_angle(b, a, c) || fortunate_angle(c, a, b);
}

i64 fortunate_sum_bruteforce(int perimeter_limit) {
    i64 answer = 0;
    for (int a = 1; a <= perimeter_limit / 3; ++a) {
        for (int b = a; a + b <= perimeter_limit; ++b) {
            const int cmax = std::min(perimeter_limit - a - b, a + b - 1);
            for (int c = b; c <= cmax; ++c) {
                if (is_fortunate_triangle(a, b, c)) {
                    answer += static_cast<i64>(a) + b + c;
                }
            }
        }
    }
    return answer;
}

void validate() {
    assert(fortunate_sum_fast(10) == 24);
    assert(fortunate_sum_fast(100) == 3331);

    for (int perimeter_limit = 5; perimeter_limit <= 220; ++perimeter_limit) {
        assert(fortunate_sum_fast(perimeter_limit) == fortunate_sum_bruteforce(perimeter_limit));
    }
}

}  // namespace

int main() {
    validate();
    std::cout << fortunate_sum_fast(10'000'000) << '\n';
    return 0;
}
