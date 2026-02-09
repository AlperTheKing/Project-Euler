#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    int perimeter_limit = 100'000;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    bool run_checkpoints = true;
};

struct PointI {
    i64 x = 0;
    i64 y = 0;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    if (options.threads <= 0) {
        options.threads = 1;
    }
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--perimeter=", options.perimeter_limit) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.perimeter_limit >= 1 && options.threads >= 1;
}

bool lex_leq(const PointI& a, const PointI& b) {
    if (a.x != b.x) {
        return a.x < b.x;
    }
    return a.y <= b.y;
}

bool point_equal(const PointI& a, const PointI& b) {
    return a.x == b.x && a.y == b.y;
}

i64 isqrt_i64(i64 x) {
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

long double dist(const PointI& a, const PointI& b) {
    const long double dx = static_cast<long double>(a.x - b.x);
    const long double dy = static_cast<long double>(a.y - b.y);
    return std::sqrt(dx * dx + dy * dy);
}

long double solve(const int perimeter_limit, int threads) {
    // From AB^2 = 4R^2 - |A+B|^2 with A+B = H-C and |H|=5,
    // each side has a lower bound sqrt(3R^2 - 10R - 25).
    // So p >= 3*sqrt(3R^2 - 10R - 25).
    const long double p = static_cast<long double>(perimeter_limit);
    const long double r_bound =
        (10.0L + std::sqrt(100.0L + 12.0L * (25.0L + (p * p) / 9.0L))) / 6.0L;
    const int r_max = static_cast<int>(std::floor(r_bound));

    threads = std::max(1, std::min(threads, r_max + 2));

    std::atomic<int> next_x{-r_max};
    constexpr int kChunk = 32;

    std::vector<long double> partial_sum(static_cast<std::size_t>(threads), 0.0L);
    std::vector<u64> partial_count(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            long double local_sum = 0.0L;
            u64 local_count = 0;

            while (true) {
                const int x_start = next_x.fetch_add(kChunk, std::memory_order_relaxed);
                if (x_start > 1) {
                    break;
                }
                const int x_end = std::min(1, x_start + kChunk - 1);

                for (int cx = x_start; cx <= x_end; ++cx) {
                    const i64 cx2 = static_cast<i64>(cx) * static_cast<i64>(cx);
                    const i64 cy_max = isqrt_i64(static_cast<i64>(r_max) * r_max - cx2);

                    for (i64 cy = -cy_max; cy <= cy_max; ++cy) {
                        const i64 r2 = cx2 + cy * cy;
                        if (r2 == 0) {
                            continue;
                        }

                        const i64 sx = 5 - static_cast<i64>(cx);
                        const i64 sy = -cy;
                        const i64 den = sx * sx + sy * sy;
                        if (den == 0) {
                            continue;
                        }

                        const i64 d = 4 * r2 - den;
                        if (d <= 0) {
                            continue;
                        }

                        const i64 g = std::gcd(std::llabs(sx), std::llabs(sy));
                        const i64 num = d * g * g;
                        if (num % den != 0) {
                            continue;
                        }

                        const i64 t2 = num / den;
                        const i64 tt = isqrt_i64(t2);
                        if (tt * tt != t2) {
                            continue;
                        }

                        const i64 px = -sy / g;
                        const i64 py = sx / g;
                        const i64 ux = tt * px;
                        const i64 uy = tt * py;

                        if (((sx + ux) & 1LL) != 0LL || ((sy + uy) & 1LL) != 0LL) {
                            continue;
                        }

                        const PointI a{(sx + ux) / 2, (sy + uy) / 2};
                        const PointI b{(sx - ux) / 2, (sy - uy) / 2};
                        const PointI c{cx, cy};

                        if (point_equal(a, b) || point_equal(a, c) || point_equal(b, c)) {
                            continue;
                        }

                        if (!lex_leq(c, a) || !lex_leq(c, b)) {
                            continue;
                        }

                        if (a.x * a.x + a.y * a.y != r2 || b.x * b.x + b.y * b.y != r2) {
                            continue;
                        }

                        if (a.x + b.x + c.x != 5 || a.y + b.y + c.y != 0) {
                            continue;
                        }

                        const i64 area2 = std::llabs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
                        if (area2 == 0) {
                            continue;
                        }

                        const long double perimeter = dist(a, b) + dist(a, c) + dist(b, c);
                        if (perimeter > static_cast<long double>(perimeter_limit) + 1e-12L) {
                            continue;
                        }

                        ++local_count;
                        local_sum += perimeter;
                    }
                }
            }

            partial_sum[static_cast<std::size_t>(t)] = local_sum;
            partial_count[static_cast<std::size_t>(t)] = local_count;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    long double total = 0.0L;
    u64 count = 0;
    for (int t = 0; t < threads; ++t) {
        total += partial_sum[static_cast<std::size_t>(t)];
        count += partial_count[static_cast<std::size_t>(t)];
    }

    (void)count;
    return total;
}

long double round4(long double x) {
    return std::round(x * 10000.0L) / 10000.0L;
}

bool run_checkpoints() {
    const long double sum50 = round4(solve(50, 1));
    if (std::fabsl(sum50 - 291.0089L) > 1e-9L) {
        std::cerr << "Checkpoint failed for perimeter=50: got " << std::setprecision(10)
                  << static_cast<double>(sum50) << '\n';
        return false;
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) {
        hw = 2;
    }
    const int multi_threads = static_cast<int>(std::min<unsigned>(hw, 8));
    const long double s1 = round4(solve(300, 1));
    const long double sm = round4(solve(300, multi_threads));
    if (std::fabsl(s1 - sm) > 1e-9L) {
        std::cerr << "Thread consistency failed for perimeter=300" << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const long double answer = round4(solve(options.perimeter_limit, options.threads));
    std::cout << std::fixed << std::setprecision(4) << static_cast<double>(answer) << '\n';
    return 0;
}
