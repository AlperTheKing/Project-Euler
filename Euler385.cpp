#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr i64 kDefaultN = 1'000'000'000LL;

struct Options {
    i64 n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct Point {
    i64 x = 0;
    i64 y = 0;

    bool operator<(const Point& other) const {
        if (x != other.x) {
            return x < other.x;
        }
        return y < other.y;
    }

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

struct TriangleKey {
    std::array<Point, 3> vertices{};

    bool operator==(const TriangleKey& other) const {
        return vertices == other.vertices;
    }
};

struct TriangleKeyHash {
    std::size_t operator()(const TriangleKey& key) const {
        std::size_t h = 0ULL;
        for (const Point& p : key.vertices) {
            const std::size_t hx = std::hash<i64>{}(p.x);
            const std::size_t hy = std::hash<i64>{}(p.y);
            h ^= hx + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
            h ^= hy + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
        }
        return h;
    }
};

struct Direction {
    int q0 = 0;
    int u = 0;
    int v = 0;
    int m = 3;
    int pell_n = 0;  // k^2 - 3 d^2 = pell_n.
    i64 d_max = 0;
    i64 k_max = 0;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            if (parsed_u64 > static_cast<u64>(std::numeric_limits<i64>::max())) {
                std::cerr << "--n is too large for this implementation.\n";
                return false;
            }
            options.n = static_cast<i64>(parsed_u64);
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    constexpr std::size_t kMinUnitsForParallel = 8ULL;

    if (!allow_multithreading || workload_units < 2ULL || workload_units < kMinUnitsForParallel) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

i64 abs_i64(const i64 x) {
    return x >= 0 ? x : -x;
}

u64 abs_to_u64(const i64 x) {
    return static_cast<u64>(x >= 0 ? x : -x);
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }

    std::string digits;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(digits.begin(), digits.end());
    return digits;
}

std::vector<std::pair<i64, i64>> generate_pell_pairs(const int pell_n,
                                                      const i64 d_max,
                                                      const i64 k_max) {
    // Problem-specific seeds for k^2 - 3 d^2 = -39 and = -3.
    std::vector<std::pair<i64, i64>> seeds;
    if (pell_n == -39) {
        seeds = {{3, 4}, {6, 5}};
    } else if (pell_n == -3) {
        seeds = {{3, 2}};
    } else {
        return {};
    }

    std::vector<std::pair<i64, i64>> all;
    for (const auto [k0, d0] : seeds) {
        i64 k = k0;
        i64 d = d0;

        while (d <= d_max && k <= k_max) {
            all.emplace_back(d, k);

            if (k > (std::numeric_limits<i64>::max() - 3LL * d) / 2LL ||
                d > (std::numeric_limits<i64>::max() - k) / 2LL) {
                break;
            }

            const i64 next_k = 2LL * k + 3LL * d;
            const i64 next_d = k + 2LL * d;
            k = next_k;
            d = next_d;
        }
    }

    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    return all;
}

std::vector<Direction> build_directions(const i64 n) {
    std::vector<Direction> directions;

    for (const int q0 : {3, 39}) {
        const int pell_n = -117 / q0;  // Derived from q0(3 d^2 - k^2) = 117.
        const int limit = static_cast<int>(std::sqrt(static_cast<long double>(q0))) + 2;

        for (int u = -limit; u <= limit; ++u) {
            for (int v = -limit; v <= limit; ++v) {
                if (u == 0 && v == 0) {
                    continue;
                }
                if (std::gcd(abs_i64(u), abs_i64(v)) != 1) {
                    continue;
                }
                if (u * u + u * v + v * v != q0) {
                    continue;
                }

                const int a = 2 * u + v;
                const int b = u + 2 * v;
                const int m = static_cast<int>(std::gcd(abs_i64(a), abs_i64(b)));
                if (m != 3) {
                    continue;
                }

                const i64 x_scale = std::max({abs_i64(u), abs_i64(v), abs_i64(u + v)});
                const i64 y_scale_raw = std::max({abs_i64(a), abs_i64(b), abs_i64(u - v)});
                if (x_scale == 0 || y_scale_raw % m != 0) {
                    continue;
                }
                const i64 y_scale = y_scale_raw / m;
                if (y_scale == 0) {
                    continue;
                }

                Direction dir;
                dir.q0 = q0;
                dir.u = u;
                dir.v = v;
                dir.m = m;
                dir.pell_n = pell_n;
                dir.d_max = n / x_scale;
                dir.k_max = n / y_scale;
                if (dir.d_max > 0 && dir.k_max > 0) {
                    directions.push_back(dir);
                }
            }
        }
    }

    std::sort(directions.begin(),
              directions.end(),
              [](const Direction& lhs, const Direction& rhs) {
                  if (lhs.q0 != rhs.q0) {
                      return lhs.q0 < rhs.q0;
                  }
                  if (lhs.u != rhs.u) {
                      return lhs.u < rhs.u;
                  }
                  return lhs.v < rhs.v;
              });
    return directions;
}

bool build_triangle_key(const Direction& dir,
                        const i64 d,
                        const i64 k,
                        const i64 n,
                        TriangleKey& out_key,
                        u64& out_area) {
    const i64 u = dir.u;
    const i64 v = dir.v;
    const i64 m = dir.m;

    const i64 a = 2LL * u + v;
    const i64 b = u + 2LL * v;
    if (a % m != 0LL || b % m != 0LL || (u - v) % m != 0LL) {
        return false;
    }

    const i64 x1 = d * u;
    const i64 x2 = d * v;
    const i64 x3 = -d * (u + v);

    const i64 y1 = (b / m) * k;
    const i64 y2 = -(a / m) * k;
    const i64 y3 = ((u - v) / m) * k;

    if (abs_i64(x1) > n || abs_i64(x2) > n || abs_i64(x3) > n || abs_i64(y1) > n ||
        abs_i64(y2) > n || abs_i64(y3) > n) {
        return false;
    }

    std::array<Point, 3> points = {{{x1, y1}, {x2, y2}, {x3, y3}}};
    std::sort(points.begin(), points.end());

    if (points[0] == points[1] || points[1] == points[2]) {
        return false;
    }

    out_key.vertices = points;
    out_area = static_cast<u64>(dir.q0) * abs_to_u64(d * k);
    return true;
}

u128 sum_triangle_areas_from_map(const std::unordered_map<TriangleKey, u64, TriangleKeyHash>& map) {
    u128 total = 0U;
    for (const auto& entry : map) {
        total += static_cast<u128>(entry.second);
    }
    return total;
}

u128 solve_fast(const i64 n,
                const bool allow_multithreading,
                const unsigned requested_threads,
                std::size_t* out_triangle_count = nullptr) {
    const std::vector<Direction> directions = build_directions(n);
    if (directions.empty()) {
        if (out_triangle_count != nullptr) {
            *out_triangle_count = 0ULL;
        }
        return 0U;
    }

    i64 max_d_n39 = 0;
    i64 max_k_n39 = 0;
    i64 max_d_n3 = 0;
    i64 max_k_n3 = 0;

    for (const Direction& dir : directions) {
        if (dir.pell_n == -39) {
            max_d_n39 = std::max(max_d_n39, dir.d_max);
            max_k_n39 = std::max(max_k_n39, dir.k_max);
        } else if (dir.pell_n == -3) {
            max_d_n3 = std::max(max_d_n3, dir.d_max);
            max_k_n3 = std::max(max_k_n3, dir.k_max);
        }
    }

    const std::vector<std::pair<i64, i64>> pairs_n39 = generate_pell_pairs(-39, max_d_n39, max_k_n39);
    const std::vector<std::pair<i64, i64>> pairs_n3 = generate_pell_pairs(-3, max_d_n3, max_k_n3);

    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, directions.size());
    std::vector<std::unordered_map<TriangleKey, u64, TriangleKeyHash>> partial_maps(
        static_cast<std::size_t>(threads));
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (unsigned t = 0U; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            auto& local = partial_maps[static_cast<std::size_t>(t)];

            for (std::size_t i = static_cast<std::size_t>(t); i < directions.size();
                 i += static_cast<std::size_t>(threads)) {
                const Direction& dir = directions[i];
                const auto& pairs = (dir.pell_n == -39) ? pairs_n39 : pairs_n3;

                for (const auto [d, k_abs] : pairs) {
                    if (d > dir.d_max || k_abs > dir.k_max) {
                        continue;
                    }

                    for (const i64 sign : {-1LL, 1LL}) {
                        TriangleKey key;
                        u64 area = 0ULL;
                        if (!build_triangle_key(dir, d, sign * k_abs, n, key, area)) {
                            continue;
                        }

                        const auto it = local.find(key);
                        if (it == local.end()) {
                            local.emplace(key, area);
                        } else if (it->second != area) {
                            std::cerr << "Internal inconsistency: same triangle got two areas.\n";
                            std::abort();
                        }
                    }
                }
            }
        });
    }

    for (std::thread& th : pool) {
        th.join();
    }

    std::unordered_map<TriangleKey, u64, TriangleKeyHash> merged;
    for (auto& local : partial_maps) {
        for (auto& kv : local) {
            const auto it = merged.find(kv.first);
            if (it == merged.end()) {
                merged.emplace(std::move(kv));
            } else if (it->second != kv.second) {
                std::cerr << "Internal inconsistency during merge.\n";
                std::abort();
            }
        }
    }

    if (out_triangle_count != nullptr) {
        *out_triangle_count = merged.size();
    }
    return sum_triangle_areas_from_map(merged);
}

u128 solve_bruteforce_small(const int n) {
    std::unordered_map<TriangleKey, u64, TriangleKeyHash> triangles;

    for (int x1 = -n; x1 <= n; ++x1) {
        for (int y1 = -n; y1 <= n; ++y1) {
            for (int x2 = -n; x2 <= n; ++x2) {
                const i64 coeff = static_cast<i64>(x1) + 2LL * static_cast<i64>(x2);
                const i64 rhs = -(2LL * static_cast<i64>(x1) + static_cast<i64>(x2)) *
                                static_cast<i64>(y1);

                const auto check_candidate = [&](const i64 y2) {
                    if (y2 < -n || y2 > n) {
                        return;
                    }

                    const i64 re = static_cast<i64>(x1) * static_cast<i64>(x1) -
                                   static_cast<i64>(y1) * static_cast<i64>(y1) +
                                   static_cast<i64>(x1) * static_cast<i64>(x2) -
                                   static_cast<i64>(y1) * y2 +
                                   static_cast<i64>(x2) * static_cast<i64>(x2) - y2 * y2;
                    if (re != 39LL) {
                        return;
                    }

                    const i64 x3 = -static_cast<i64>(x1) - static_cast<i64>(x2);
                    const i64 y3 = -static_cast<i64>(y1) - y2;
                    if (abs_i64(x3) > n || abs_i64(y3) > n) {
                        return;
                    }

                    std::array<Point, 3> points = {{{x1, y1}, {x2, y2}, {x3, y3}}};
                    std::sort(points.begin(), points.end());
                    if (points[0] == points[1] || points[1] == points[2]) {
                        return;
                    }

                    const i64 area2 = abs_i64((points[1].x - points[0].x) * (points[2].y - points[0].y) -
                                              (points[1].y - points[0].y) * (points[2].x - points[0].x));
                    if (area2 == 0 || (area2 & 1LL) != 0LL) {
                        return;
                    }

                    TriangleKey key;
                    key.vertices = points;
                    triangles[key] = static_cast<u64>(area2 / 2LL);
                };

                if (coeff == 0LL) {
                    if (rhs != 0LL) {
                        continue;
                    }
                    for (int y2 = -n; y2 <= n; ++y2) {
                        check_candidate(static_cast<i64>(y2));
                    }
                } else {
                    if (rhs % coeff != 0LL) {
                        continue;
                    }
                    check_candidate(rhs / coeff);
                }
            }
        }
    }

    return sum_triangle_areas_from_map(triangles);
}

bool run_checkpoints(const bool allow_multithreading, const unsigned requested_threads) {
    struct Checkpoint {
        int n;
        u64 expected;
    };

    const std::vector<Checkpoint> known = {
        {8, 72ULL},
        {10, 252ULL},
        {100, 34'632ULL},
        {1'000, 3'529'008ULL},
    };

    for (const Checkpoint& cp : known) {
        const u128 got = solve_fast(cp.n, allow_multithreading, requested_threads);
        if (got != static_cast<u128>(cp.expected)) {
            std::cerr << "Checkpoint failed for A(" << cp.n << "): expected " << cp.expected
                      << ", got " << to_string_u128(got) << '\n';
            return false;
        }
    }

    const u128 brute20 = solve_bruteforce_small(20);
    const u128 fast20 = solve_fast(20, false, 1U);
    if (brute20 != fast20) {
        std::cerr << "Bruteforce cross-check failed for n=20: brute=" << to_string_u128(brute20)
                  << ", fast=" << to_string_u128(fast20) << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.n < 0) {
        std::cerr << "--n must be non-negative.\n";
        return 1;
    }

    if (options.run_checkpoints &&
        !run_checkpoints(options.allow_multithreading, options.requested_threads)) {
        return 1;
    }

    std::size_t triangle_count = 0ULL;
    const u128 answer = solve_fast(options.n,
                                   options.allow_multithreading,
                                   options.requested_threads,
                                   &triangle_count);

    std::cout << to_string_u128(answer) << '\n';
    std::cerr << "Triangles counted: " << triangle_count << '\n';
    return 0;
}
