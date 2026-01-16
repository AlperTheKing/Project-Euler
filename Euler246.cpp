#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = long long;
using u64 = unsigned long long;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr i64 kMx = -2000;
constexpr i64 kMy = 1500;
constexpr i64 kGx = 8000;
constexpr i64 kGy = 1500;
constexpr i64 kRadius = 15000;

u64 isqrt_u128(u128 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while (static_cast<u128>(r + 1) * (r + 1) <= x) ++r;
    while (static_cast<u128>(r) * r > x) --r;
    return r;
}

struct EllipseData {
    i64 a2 = 0;
    i64 b2 = 0;
    i64 r2_circle = 0;
    i128 ab2 = 0;
};

i64 min_y_outside(i64 x, const EllipseData& e) {
    u128 lhs = static_cast<u128>(e.b2) * x * x;
    u128 rhs = static_cast<u128>(e.a2) * e.b2;
    if (lhs > rhs) return 0;
    if (lhs == rhs) return 1;
    u128 rem = rhs - lhs;
    u128 q = rem / static_cast<u128>(e.a2);
    i64 y = static_cast<i64>(isqrt_u128(q));
    while (static_cast<u128>(e.a2) * y * y <= rem) ++y;
    return y;
}

i128 angle_ineq(i64 x, i128 u, const EllipseData& e) {
    // Positive when the acute angle between tangents is greater than 45 degrees.
    i128 x2 = static_cast<i128>(x) * x;
    i128 term1 = 4 * (static_cast<i128>(e.b2) * x2 +
                      static_cast<i128>(e.a2) * u - e.ab2);
    i128 t = static_cast<i128>(e.a2) + e.b2 - x2 - u;
    return term1 - t * t;
}

i64 max_y_outside(i64 x, const EllipseData& e) {
    i128 x2 = static_cast<i128>(x) * x;
    i128 k = static_cast<i128>(e.a2) + e.b2 - x2;
    i128 b = 4 * static_cast<i128>(e.a2) + 2 * k;
    i128 c = 4 * (static_cast<i128>(e.b2) * x2 - e.ab2) - k * k;
    i128 d = b * b + 4 * c;
    if (d <= 0) return -1;
    u128 sqrt_d = isqrt_u128(static_cast<u128>(d));
    i128 u_max = (b + static_cast<i128>(sqrt_d)) / 2;
    if (u_max < 0) return -1;
    while (u_max >= 0 && angle_ineq(x, u_max, e) <= 0) --u_max;
    if (u_max < 0) return -1;
    i64 y = static_cast<i64>(isqrt_u128(static_cast<u128>(u_max)));
    while (y >= 0 && angle_ineq(x, static_cast<i128>(y) * y, e) <= 0) --y;
    return y;
}

bool qualifies(i64 x, i64 y, const EllipseData& e) {
    i128 x2 = static_cast<i128>(x) * x;
    i128 y2 = static_cast<i128>(y) * y;
    if (static_cast<i128>(e.b2) * x2 + static_cast<i128>(e.a2) * y2 <= e.ab2) {
        return false;
    }
    i128 r2 = x2 + y2;
    if (r2 <= e.r2_circle) return true;
    return angle_ineq(x, y2, e) > 0;
}

i64 count_quadrant(const EllipseData& e, int threads) {
    i64 x_limit = static_cast<i64>(isqrt_u128(static_cast<u128>(e.a2) +
                                             static_cast<u128>(6) * e.b2)) +
                  2;
    int use_threads = std::max(1, threads);
    use_threads = std::min<i64>(use_threads, x_limit + 1);
    std::vector<i64> locals(use_threads, 0);

    auto worker = [&](int tid, i64 start, i64 end) {
        i64 local = 0;
        for (i64 x = start; x < end; ++x) {
            i64 y_min = min_y_outside(x, e);
            i64 y_circle = -1;
            i64 x2 = x * x;
            if (x2 <= e.r2_circle) {
                // Inside the director circle (r^2 <= a^2 + b^2), the tangent angle is obtuse.
                y_circle = static_cast<i64>(isqrt_u128(e.r2_circle - x2));
            }
            if (y_circle >= y_min) {
                local += y_circle - y_min + 1;
            }
            i64 y_start = std::max(y_min, y_circle + 1);
            i64 y_max = max_y_outside(x, e);
            if (y_max >= y_start) {
                local += y_max - y_start + 1;
            }
        }
        locals[tid] = local;
    };

    i64 total_x = x_limit + 1;
    i64 chunk = (total_x + use_threads - 1) / use_threads;
    std::vector<std::thread> pool;
    pool.reserve(use_threads);
    for (int t = 0; t < use_threads; ++t) {
        i64 start = t * chunk;
        i64 end = std::min(start + chunk, total_x);
        if (start >= end) break;
        pool.emplace_back(worker, t, start, end);
    }
    for (auto& th : pool) th.join();

    i64 sum = 0;
    for (i64 v : locals) sum += v;
    return sum;
}

i64 count_points(const EllipseData& e, int threads) {
    i64 count_q1 = count_quadrant(e, threads);

    i64 x_limit = static_cast<i64>(isqrt_u128(static_cast<u128>(e.a2) +
                                             static_cast<u128>(6) * e.b2)) +
                  2;
    i64 y_limit = static_cast<i64>(isqrt_u128(static_cast<u128>(e.b2) +
                                             static_cast<u128>(6) * e.a2)) +
                  2;

    i64 count_x_axis = 0;
    for (i64 x = 0; x <= x_limit; ++x) {
        if (qualifies(x, 0, e)) ++count_x_axis;
    }
    i64 count_y_axis = 0;
    for (i64 y = 0; y <= y_limit; ++y) {
        if (qualifies(0, y, e)) ++count_y_axis;
    }
    i64 origin = qualifies(0, 0, e) ? 1 : 0;

    return 4 * count_q1 - 2 * count_x_axis - 2 * count_y_axis + origin;
}

EllipseData build_ellipse(i64 mx, i64 my, i64 gx, i64 gy, i64 r) {
    if (my != gy) {
        std::cerr << "Ellipse axis is not horizontal.\n";
        std::exit(1);
    }
    if ((mx + gx) % 2 != 0 || (my + gy) % 2 != 0) {
        std::cerr << "Ellipse center is not on lattice.\n";
        std::exit(1);
    }
    i64 d = std::llabs(gx - mx);
    if (r <= d) {
        std::cerr << "Radius must exceed focal distance.\n";
        std::exit(1);
    }
    if ((r % 2) != 0 || (d % 2) != 0) {
        std::cerr << "Radius and focal distance must be even.\n";
        std::exit(1);
    }
    i64 a = r / 2;
    i64 c = d / 2;
    i64 a2 = a * a;
    i64 b2 = a2 - c * c;
    if (b2 <= 0) {
        std::cerr << "Invalid ellipse parameters.\n";
        std::exit(1);
    }
    EllipseData e;
    e.a2 = a2;
    e.b2 = b2;
    e.r2_circle = a2 + b2;
    e.ab2 = static_cast<i128>(a2) * b2;
    return e;
}

void run_validation() {
    struct Check {
        i64 r;
        i64 d;
        i64 expected;
    };
    const Check checks[] = {
        {6, 2, 154},
        {10, 4, 418},
        {14, 6, 810},
        {18, 8, 1334},
    };
    for (const auto& chk : checks) {
        i64 a = chk.r / 2;
        i64 c = chk.d / 2;
        EllipseData e;
        e.a2 = a * a;
        e.b2 = e.a2 - c * c;
        e.r2_circle = e.a2 + e.b2;
        e.ab2 = static_cast<i128>(e.a2) * e.b2;
        i64 got = count_points(e, 1);
        if (got != chk.expected) {
            std::cerr << "Validation failed for r=" << chk.r << " d=" << chk.d
                      << ": got " << got << ", expected " << chk.expected << '\n';
            std::exit(1);
        }
    }
}

void usage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " [-t THREADS] [--no-validate]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-t" && i + 1 < argc) {
            threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--validate") {
            validate = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (validate) run_validation();

    EllipseData ellipse = build_ellipse(kMx, kMy, kGx, kGy, kRadius);
    i64 answer = count_points(ellipse, threads);
    std::cout << answer << '\n';
    return 0;
}
