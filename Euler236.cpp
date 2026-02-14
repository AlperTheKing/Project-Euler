#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;

struct Fraction {
    i64 u = 0;
    i64 v = 1;
};

struct Options {
    int threads = 0;
    bool run_checkpoints = true;
};

constexpr int NUM_PRODUCTS = 5;
constexpr std::array<i64, NUM_PRODUCTS> A = {5248, 1312, 2624, 5760, 3936};
constexpr std::array<i64, NUM_PRODUCTS> B = {640, 1888, 3776, 3776, 5664};
constexpr int GROUPS = 3;
constexpr std::array<std::array<int, 3>, GROUPS> GROUP_PRODUCTS = {{{0, -1, -1}, {1, 2, 4}, {3, -1, -1}}};
constexpr std::array<int, GROUPS> GROUP_SIZE = {1, 3, 1};
constexpr std::array<i64, GROUPS> GROUP_B_NUM = {5, 59, 59};
constexpr std::array<i64, GROUPS> GROUP_A_DEN = {41, 41, 90};
constexpr i64 SA = 18880;
constexpr i64 SB = 15744;

i64 abs_i64(const i64 x) {
    return x < 0 ? -x : x;
}

Fraction reduce_fraction(i64 u, i64 v) {
    if (v < 0) {
        u = -u;
        v = -v;
    }
    const i64 g = std::gcd(abs_i64(u), abs_i64(v));
    return {u / g, v / g};
}

bool fraction_equal(const Fraction& a, const Fraction& b) {
    return a.u == b.u && a.v == b.v;
}

bool fraction_less(const Fraction& a, const Fraction& b) {
    return static_cast<i128>(a.u) * static_cast<i128>(b.v) < static_cast<i128>(b.u) * static_cast<i128>(a.v);
}

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
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.threads >= 0;
}

i64 floor_div(i64 a, i64 b) {
    i64 q = a / b;
    i64 r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) {
        --q;
    }
    return q;
}

i64 ceil_div(i64 a, i64 b) {
    return -floor_div(-a, b);
}

i64 extended_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = (a >= 0) ? 1 : -1;
        y = 0;
        return abs_i64(a);
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

bool tighten_range(i64 p, i64 q, i64 low, i64 high, i64& t_low, i64& t_high) {
    if (low > high) {
        return false;
    }
    if (q == 0) {
        return (p >= low && p <= high);
    }

    i64 l = 0;
    i64 r = 0;
    if (q > 0) {
        l = ceil_div(low - p, q);
        r = floor_div(high - p, q);
    } else {
        l = ceil_div(high - p, q);
        r = floor_div(low - p, q);
    }
    if (l > r) {
        return false;
    }
    if (l > t_low) {
        t_low = l;
    }
    if (r < t_high) {
        t_high = r;
    }
    return t_low <= t_high;
}

bool exists_bounded_linear(i64 c1, i64 l1, i64 h1, i64 c2, i64 l2, i64 h2, i64 rhs) {
    if (l1 > h1 || l2 > h2) {
        return false;
    }
    if (c1 == 0 && c2 == 0) {
        return rhs == 0;
    }
    if (c1 == 0) {
        if (rhs % c2 != 0) {
            return false;
        }
        const i64 s2 = rhs / c2;
        return s2 >= l2 && s2 <= h2;
    }
    if (c2 == 0) {
        if (rhs % c1 != 0) {
            return false;
        }
        const i64 s1 = rhs / c1;
        return s1 >= l1 && s1 <= h1;
    }

    i64 x0 = 0;
    i64 y0 = 0;
    const i64 g = extended_gcd(c1, c2, x0, y0);
    if (rhs % g != 0) {
        return false;
    }

    const i64 scale = rhs / g;
    const i64 step1 = c2 / g;
    const i64 step2 = -c1 / g;

    const i128 base1_128 = static_cast<i128>(x0) * static_cast<i128>(scale);
    const i128 base2_128 = static_cast<i128>(y0) * static_cast<i128>(scale);
    if (base1_128 < static_cast<i128>(std::numeric_limits<i64>::min()) ||
        base1_128 > static_cast<i128>(std::numeric_limits<i64>::max()) ||
        base2_128 < static_cast<i128>(std::numeric_limits<i64>::min()) ||
        base2_128 > static_cast<i128>(std::numeric_limits<i64>::max())) {
        return false;
    }
    const i64 base1 = static_cast<i64>(base1_128);
    const i64 base2 = static_cast<i64>(base2_128);

    i64 t_low = std::numeric_limits<i64>::min() / 4;
    i64 t_high = std::numeric_limits<i64>::max() / 4;
    if (!tighten_range(base1, step1, l1, h1, t_low, t_high)) {
        return false;
    }
    if (!tighten_range(base2, step2, l2, h2, t_low, t_high)) {
        return false;
    }
    return t_low <= t_high;
}

std::vector<Fraction> generate_candidates() {
    std::vector<Fraction> out;
    out.reserve(static_cast<std::size_t>(A[0] * B[0]));

    for (i64 a = 1; a <= A[0]; ++a) {
        for (i64 b = 1; b <= B[0]; ++b) {
            Fraction m = reduce_fraction(41 * b, 5 * a);
            if (m.u > m.v) {
                out.push_back(m);
            }
        }
    }

    std::sort(out.begin(), out.end(), fraction_less);
    out.erase(std::unique(out.begin(), out.end(), fraction_equal), out.end());
    return out;
}

bool candidate_works(const Fraction& m) {
    if (m.u <= m.v) {
        return false;
    }

    std::array<i64, GROUPS> min_sum = {1, 3, 1};
    std::array<i64, GROUPS> max_sum = {0, 0, 0};
    std::array<i64, GROUPS> coef = {0, 0, 0};

    for (int g = 0; g < GROUPS; ++g) {
        const Fraction group_ratio = reduce_fraction(m.u * GROUP_B_NUM[g], m.v * GROUP_A_DEN[g]);
        const i64 N = group_ratio.u;
        const i64 D = group_ratio.v;

        i64 upper_sum = 0;
        for (int j = 0; j < GROUP_SIZE[g]; ++j) {
            const int idx = GROUP_PRODUCTS[g][j];
            const i64 upper = std::min(A[idx] / D, B[idx] / N);
            if (upper < 1) {
                return false;
            }
            upper_sum += upper;
        }
        max_sum[g] = upper_sum;
        coef[g] = m.v * SB * D - m.u * SA * N;
    }

    for (i64 s0 = min_sum[0]; s0 <= max_sum[0]; ++s0) {
        const i64 rhs = -(coef[0] * s0);
        if (exists_bounded_linear(coef[1], min_sum[1], max_sum[1], coef[2], min_sum[2], max_sum[2], rhs)) {
            return true;
        }
    }

    return false;
}

std::vector<Fraction> solve_all(const std::vector<Fraction>& candidates, int threads) {
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
        if (threads <= 0) {
            threads = 1;
        }
    }

    const int n = static_cast<int>(candidates.size());
    threads = std::max(1, std::min(threads, n));

    std::vector<std::vector<Fraction>> local(static_cast<std::size_t>(threads));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            for (int i = t; i < n; i += threads) {
                if (candidate_works(candidates[static_cast<std::size_t>(i)])) {
                    local[static_cast<std::size_t>(t)].push_back(candidates[static_cast<std::size_t>(i)]);
                }
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }

    std::vector<Fraction> solutions;
    for (auto& chunk : local) {
        solutions.insert(solutions.end(), chunk.begin(), chunk.end());
    }
    std::sort(solutions.begin(), solutions.end(), fraction_less);
    solutions.erase(std::unique(solutions.begin(), solutions.end(), fraction_equal), solutions.end());
    return solutions;
}

bool run_checkpoints(const std::vector<Fraction>& solutions) {
    if (solutions.size() != 35U) {
        std::cerr << "Checkpoint failed for number of valid m values" << '\n';
        return false;
    }
    if (!fraction_equal(solutions.front(), Fraction{1476, 1475})) {
        std::cerr << "Checkpoint failed for smallest valid m value" << '\n';
        return false;
    }
    if (!candidate_works(Fraction{1476, 1475})) {
        std::cerr << "Checkpoint failed for explicit smallest-m verification" << '\n';
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

    const std::vector<Fraction> candidates = generate_candidates();
    const std::vector<Fraction> solutions = solve_all(candidates, options.threads);

    if (solutions.empty()) {
        std::cerr << "No solution found" << '\n';
        return 2;
    }
    if (options.run_checkpoints && !run_checkpoints(solutions)) {
        return 3;
    }

    const Fraction best = solutions.back();
    std::cout << best.u << '/' << best.v << '\n';
    return 0;
}
