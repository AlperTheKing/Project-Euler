#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;
using u64 = std::uint64_t;

template <typename T>
bool mod_inv_in_range(T a, T m, T& inv) {
    a %= m;
    if (a < 0) a += m;
    T x = a;
    T y = m;
    T vx = 1;
    T vy = 0;
    while (x) {
        T k = y / x;
        y %= x;
        vy -= k * vx;
        std::swap(x, y);
        std::swap(vx, vy);
    }
    if (y != 1) return false;
    inv = (vy < 0) ? (vy + m) : vy;
    return true;
}

u64 solve_f(i64 n) {
    for (int steps = 0; ; ++steps) {
        std::pair<int, i64> ans{steps, n + 1};
        const int d = (steps + 1) / 2;

        auto dfs = [&](auto&& self, i64 a, i64 b, i64 ca, i64 cb, int cur_depth) -> void {
            if (a > b) {
                std::swap(a, b);
                std::swap(ca, cb);
            }

            if (cb < 0) return;
            if (ca < 0) {
                const i64 k = (-ca + b - 1) / b;
                ca += k * b;
                cb -= k * a;
                if (cb < 0) return;
            }
            if (a * ca + b * cb != n) return;

            if (cur_depth == d) {
                if ((steps & 1) != 0) {
                    a *= 2;
                    b *= 2;
                    ca *= 2;
                    cb *= 2;
                    b -= a / 2;
                    ca += cb / 2;
                    if (a * ca + b * cb != 4 * n) return;
                }

                const i128 norm = static_cast<i128>(a) * a + static_cast<i128>(b) * b;
                if (norm < n) return;

                i128 cross = static_cast<i128>(ca) * b - static_cast<i128>(cb) * a;
                {
                    const i128 k = cross / norm;
                    ca -= static_cast<i64>(k * b);
                    cb += static_cast<i64>(k * a);
                    cross -= k * norm;
                }
                if (cross < 0) {
                    ca += b;
                    cb -= a;
                    cross += norm;
                }

                auto check = [&]() {
                    if (ca < 0 || cb < 0) return;
                    if (static_cast<i128>(ca) * ca + static_cast<i128>(cb) * cb > norm) return;

                    i64 x = ca;
                    i64 y = cb;
                    i64 vx = a;
                    i64 vy = b;

                    if ((steps & 1) != 0) {
                        if ((vx & 1) != 0) return;
                        if ((y & 1) != 0) return;
                        x -= y / 2;
                        vy += vx / 2;
                        if ((vy & 1) != 0) return;
                        if ((x & 1) != 0) return;
                        x /= 2;
                        y /= 2;
                        vx /= 2;
                        vy /= 2;
                    }

                    if (x * vx + y * vy != n) return;
                    int ops = 0;
                    for (; ops <= steps - d && x != 0 && y != 0; ++ops) {
                        if (x > y) {
                            std::swap(x, y);
                            std::swap(vx, vy);
                        }
                        y -= x;
                        vx += vy;
                    }
                    if (ops > steps - d) return;

                    i64 o = vx + vy - n;
                    o %= n;
                    if (o < 0) o += n;

                    i64 inv_o = 0;
                    if (!mod_inv_in_range(o, n, inv_o)) return;

                    ans = std::min(ans, std::pair<int, i64>{cur_depth + ops,
                                                            std::min({o, n - o, inv_o, n - inv_o})});
                };

                check();
                ca -= b;
                cb += a;
                check();
                return;
            }

            {
                i64 x = a;
                i64 y = b;
                for (int z = cur_depth; z < steps / 2; ++z) {
                    if (x > y) std::swap(x, y);
                    x += y;
                    std::swap(x, y);
                }
                if (x > y) std::swap(x, y);
                if ((steps & 1) != 0) {
                    if (5 * static_cast<i128>(y) * y / 4 + static_cast<i128>(x) * y +
                            static_cast<i128>(x) * x <
                        n) {
                        return;
                    }
                } else {
                    if (static_cast<i128>(x) * x + static_cast<i128>(y) * y < n) return;
                }
            }

            self(self, b, a + b, cb - ca, ca, cur_depth + 1);
            if (0 < a && a < b) {
                self(self, a, a + b, ca - cb, cb, cur_depth + 1);
            }
        };

        dfs(dfs, 0, 1, 0, n, 0);

        if (ans.second <= n) return static_cast<u64>(ans.second);
    }
}

void run_validations() {
    assert(solve_f(7) == 2);
    assert(solve_f(89) == 34);
    assert(solve_f(8'191) == 1'856);
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve_f(1'000'000'000'039LL) << '\n';
    return 0;
}
