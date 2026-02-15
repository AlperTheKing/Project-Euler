#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>

using i64 = long long;
using i128 = __int128_t;
using ResidueTable = std::array<std::array<std::vector<int>, 304>, 8>;

static std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }
    bool neg = x < 0;
    if (neg) {
        x = -x;
    }
    std::string s;
    while (x > 0) {
        int d = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    if (neg) {
        s.push_back('-');
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static int class_index_from_g(int g) {
    switch (g) {
        case 1: return 0;
        case 4: return 1;
        case 16: return 2;
        case 19: return 3;
        case 64: return 4;
        case 76: return 5;
        case 304: return 6;
        case 1216: return 7;
        default: return -1;
    }
}

static int trailing_zeros_mod16(int x_mod16) {
    x_mod16 &= 15;
    if (x_mod16 == 0) {
        return 4;
    }
    int c = 0;
    while ((x_mod16 & 1) == 0) {
        x_mod16 >>= 1;
        ++c;
    }
    return c;
}

static ResidueTable build_residue_table() {
    ResidueTable table;

    for (int vm = 0; vm < 304; ++vm) {
        for (int um = 0; um < 304; ++um) {
            const int m16 = (5 * um - 3 * vm) & 15;
            const int n16 = (3 * um - 5 * vm) & 15;

            const int tz = std::min(std::min(trailing_zeros_mod16(m16), trailing_zeros_mod16(n16)), 3);
            const int g2 = (tz == 0 ? 1 : (tz == 1 ? 4 : (tz == 2 ? 16 : 64)));

            int m19 = (5 * um - 3 * vm) % 19;
            if (m19 < 0) {
                m19 += 19;
            }
            const bool e19 = (m19 == 0);
            const int g = g2 * (e19 ? 19 : 1);

            const int idx = class_index_from_g(g);
            assert(idx >= 0);
            table[idx][vm].push_back(um);
        }
    }

    return table;
}

static const ResidueTable& get_residue_table() {
    static const ResidueTable table = build_residue_table();
    return table;
}

static i128 accumulate_stride(i64 N, const std::array<int, 8>& vmax, int tid, int step) {
    constexpr std::array<int, 8> G = {1, 4, 16, 19, 64, 76, 304, 1216};
    const auto& residue_table = get_residue_table();

    const long double lo = 5.0L / 3.0L;
    const long double hi = (103.0L - 4.0L * std::sqrt(19.0L)) / 45.0L;

    i128 ans = 0;

    for (int ci = 0; ci < 8; ++ci) {
        const int g = G[ci];
        for (int v = tid + 1; v <= vmax[ci]; v += step) {
            const i64 umin = static_cast<i64>(std::floor(lo * static_cast<long double>(v))) + 1;
            const i64 umax = static_cast<i64>(std::floor(hi * static_cast<long double>(v)));
            if (umin > umax) {
                continue;
            }

            const auto& residues = residue_table[ci][v % 304];
            for (int r : residues) {
                i64 first = umin;
                const int delta = (r - static_cast<int>(first % 304) + 304) % 304;
                first += delta;

                for (i64 u = first; u <= umax; u += 304) {
                    if (std::gcd(u, static_cast<i64>(v)) != 1) {
                        continue;
                    }

                    const i64 uu = u;
                    const i64 vv = v;

                    const i64 z0 = 32LL * uu * uu - 176LL * uu * vv + 240LL * vv * vv;
                    if (z0 > N * static_cast<i64>(g)) {
                        continue;
                    }

                    const i64 x0 = 15LL * uu * uu - 34LL * uu * vv + 15LL * vv * vv;
                    const i64 y0 = 105LL * uu * uu - 446LL * uu * vv + 473LL * vv * vv;

                    const i64 x = x0 / g;
                    const i64 y = y0 / g;
                    const i64 z = z0 / g;

                    ans += static_cast<i128>(x) + static_cast<i128>(y) + static_cast<i128>(z);
                }
            }
        }
    }

    return ans;
}

struct ThreadTask785 {
    i64 N = 0;
    std::array<int, 8> vmax{};
    int tid = 0;
    int step = 1;
    i128 partial = 0;
};

static void* thread_worker_785(void* arg) {
    auto* t = static_cast<ThreadTask785*>(arg);
    t->partial = accumulate_stride(t->N, t->vmax, t->tid, t->step);
    return nullptr;
}

static i128 solve_fast(i64 N) {
    constexpr std::array<int, 8> G = {1, 4, 16, 19, 64, 76, 304, 1216};
    const long double hi = (103.0L - 4.0L * std::sqrt(19.0L)) / 45.0L;
    const long double cmin = 32.0L * hi * hi - 176.0L * hi + 240.0L;

    std::array<int, 8> vmax{};
    for (int i = 0; i < 8; ++i) {
        const long double lim = std::sqrt((static_cast<long double>(N) * static_cast<long double>(G[i])) / cmin);
        vmax[i] = static_cast<int>(std::floor(lim)) + 3;
    }

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    int thread_count = (cpu_count > 1 ? static_cast<int>(cpu_count) : 1);
    if (thread_count > 16) {
        thread_count = 16;
    }
    if (N < 1'000'000LL || thread_count <= 1) {
        return accumulate_stride(N, vmax, 0, 1);
    }

    std::vector<pthread_t> threads((std::size_t)thread_count);
    std::vector<ThreadTask785> tasks((std::size_t)thread_count);

    for (int t = 0; t < thread_count; ++t) {
        tasks[(std::size_t)t].N = N;
        tasks[(std::size_t)t].vmax = vmax;
        tasks[(std::size_t)t].tid = t;
        tasks[(std::size_t)t].step = thread_count;
        const int rc = ::pthread_create(&threads[(std::size_t)t], nullptr, thread_worker_785, &tasks[(std::size_t)t]);
        assert(rc == 0);
    }

    i128 ans = 0;
    for (int t = 0; t < thread_count; ++t) {
        const int rc = ::pthread_join(threads[(std::size_t)t], nullptr);
        assert(rc == 0);
        ans += tasks[(std::size_t)t].partial;
    }
    return ans;
}

static i128 solve_bruteforce(int N) {
    i128 ans = 0;
    for (int x = 1; x <= N; ++x) {
        for (int y = x; y <= N; ++y) {
            const i64 A = 15;
            const i64 B = -34LL * (x + y);
            const i64 C = 15LL * (1LL * x * x + 1LL * y * y) - 34LL * x * y;
            const i64 D = B * B - 4LL * A * C;
            const i64 t = static_cast<i64>(std::sqrt(static_cast<long double>(D)));
            if (t * t != D) {
                continue;
            }
            for (i64 sgn : {-1LL, 1LL}) {
                const i64 num = -B + sgn * t;
                if (num % (2LL * A) != 0) {
                    continue;
                }
                const i64 z = num / (2LL * A);
                if (z < y || z > N) {
                    continue;
                }
                const i64 lhs = 15LL * (1LL * x * x + 1LL * y * y + z * z);
                const i64 rhs = 34LL * (1LL * x * y + 1LL * y * z + 1LL * z * x);
                if (lhs != rhs) {
                    continue;
                }
                if (std::gcd(std::gcd(x, y), static_cast<int>(z)) != 1) {
                    continue;
                }
                ans += static_cast<i128>(x + y + z);
            }
        }
    }
    return ans;
}

int main() {
    assert(to_string_i128(solve_fast(100)) == "184");
    assert(solve_fast(200) == solve_bruteforce(200));

    const i128 ans = solve_fast(1'000'000'000LL);
    std::cout << to_string_i128(ans) << '\n';
    return 0;
}
