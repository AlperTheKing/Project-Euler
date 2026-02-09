#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    int limit = 100000;
    unsigned threads = std::thread::hardware_concurrency();
    bool run_checkpoints = true;
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

bool parse_unsigned_after_prefix(const std::string& arg, const std::string& prefix, unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    unsigned parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10U + static_cast<unsigned>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.limit < 3) {
        return false;
    }
    if (options.threads == 0) {
        options.threads = 1;
    }
    return true;
}

long double be_from_geometry(const int a, const int b, const int c) {
    // C=(0,0), A=(b,0), B determined by side lengths.
    const long double bx = (static_cast<long double>(a) * a + static_cast<long double>(b) * b -
                            static_cast<long double>(c) * c) /
                           (2.0L * b);
    const long double by2 = static_cast<long double>(a) * a - bx * bx;
    const long double by = -std::sqrt(std::max(0.0L, by2));

    // Angle bisector at C points along unit(CA)+unit(CB).
    const long double ux = 1.0L + bx / a;
    const long double uy = by / a;

    // Tangent at C to circumcircle is perpendicular to OC.
    // Circumcenter O satisfies OA=OC and OB=OC.
    // With C at origin and A on x-axis: O=(b/2, oy).
    const long double oy =
        (bx * bx + by * by - b * bx) / (2.0L * by);  // from OB^2 = OC^2.
    const long double ox = b / 2.0L;

    // Tangent direction is perpendicular to O vector.
    const long double tx = -oy;
    const long double ty = ox;

    // Solve C + s*u = B + t*tangent, then BE = |t|*|tangent|.
    const long double det = ux * (-ty) - uy * (-tx);
    const long double rhsx = bx;
    const long double rhsy = by;

    const long double t = (ux * rhsy - uy * rhsx) / det;
    const long double tangent_len = std::sqrt(tx * tx + ty * ty);
    return std::fabsl(t) * tangent_len;
}

u64 count_range(const int limit, const int a_start, const int a_end) {
    u64 total = 0;

    for (int a = a_start; a <= a_end; ++a) {
        const int b_max = (limit - a) / 2;
        for (int b = a; b <= b_max; ++b) {
            const int s = a + b;
            const int c_max = std::min(s - 1, limit - s);
            if (c_max < b) {
                continue;
            }

            const int g = std::gcd(a, b);
            const int step = s / g;
            total += static_cast<u64>(c_max / step - (b - 1) / step);
        }
    }

    return total;
}

u64 solve_fast(const int limit, unsigned threads) {
    const int a_max = limit / 3;
    if (threads <= 1 || a_max < 5000) {
        return count_range(limit, 1, a_max);
    }

    const unsigned use_threads = std::min<unsigned>(threads, static_cast<unsigned>(a_max));
    const int chunk = (a_max + static_cast<int>(use_threads) - 1) / static_cast<int>(use_threads);

    std::vector<std::thread> pool;
    std::vector<u64> partial(use_threads, 0);
    pool.reserve(use_threads);

    for (unsigned t = 0; t < use_threads; ++t) {
        const int start = static_cast<int>(t) * chunk + 1;
        const int end = std::min(a_max, start + chunk - 1);
        if (start > end) {
            continue;
        }
        pool.emplace_back([&, start, end, t]() {
            partial[t] = count_range(limit, start, end);
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u64 total = 0;
    for (u64 part : partial) {
        total += part;
    }
    return total;
}

u64 solve_bruteforce(const int limit) {
    u64 total = 0;
    for (int a = 1; a <= limit / 3; ++a) {
        for (int b = a; b <= (limit - a) / 2; ++b) {
            const int c_max = std::min(a + b - 1, limit - a - b);
            for (int c = b; c <= c_max; ++c) {
                if ((static_cast<i64>(a) * c) % (a + b) == 0) {
                    ++total;
                }
            }
        }
    }
    return total;
}

bool run_checkpoints() {
    {
        // Geometry validation for the BE formula on a few integer triangles.
        struct Tri {
            int a;
            int b;
            int c;
        };
        const std::vector<Tri> tests = {
            {3, 4, 5},
            {5, 5, 8},
            {7, 8, 9},
            {9, 12, 14},
            {11, 13, 15},
        };

        for (const auto& tri : tests) {
            const long double be_geom = be_from_geometry(tri.a, tri.b, tri.c);
            const long double be_formula =
                static_cast<long double>(tri.a) * tri.c / static_cast<long double>(tri.a + tri.b);
            if (std::fabsl(be_geom - be_formula) > 1e-9L) {
                std::cerr << "Geometry checkpoint failed for (" << tri.a << ',' << tri.b << ','
                          << tri.c << ")" << '\n';
                return false;
            }
        }
    }

    for (int limit : {150, 300, 600}) {
        const u64 fast = solve_fast(limit, 1);
        const u64 brute = solve_bruteforce(limit);
        if (fast != brute) {
            std::cerr << "Counting checkpoint failed for limit=" << limit << ": fast=" << fast
                      << ", brute=" << brute << '\n';
            return false;
        }
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

    const u64 answer = solve_fast(options.limit, options.threads);
    std::cout << answer << '\n';
    return 0;
}
