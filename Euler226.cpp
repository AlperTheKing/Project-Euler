#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    long double eps = 2e-9L;
    int min_level = 22;
    int max_level = 28;
    int threads = 0;
    bool run_checkpoints = true;
};

bool parse_ld_after_prefix(const std::string& arg, const std::string& prefix, long double& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    try {
        value = std::stold(tail);
    } catch (...) {
        return false;
    }
    return value > 0;
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
        if (parse_ld_after_prefix(arg, "--eps=", options.eps) ||
            parse_int_after_prefix(arg, "--min-level=", options.min_level) ||
            parse_int_after_prefix(arg, "--max-level=", options.max_level) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.eps > 0.0L && options.min_level >= 3 && options.max_level >= options.min_level &&
           options.max_level <= 30 && options.threads >= 0;
}

u64 prefix_popcount(u64 n) {
    u64 total = 0;
    for (int bit = 0; bit < 63; ++bit) {
        const u64 half = 1ULL << bit;
        if (half >= n) {
            break;
        }
        const u64 cycle = half << 1U;
        const u64 full = n / cycle;
        const u64 rem = n % cycle;
        total += full * half;
        if (rem > half) {
            total += rem - half;
        }
    }
    return total;
}

long double takagi_dyadic(const u64 k, const int level) {
    const long double numer = static_cast<long double>(level) * static_cast<long double>(k) -
                              2.0L * static_cast<long double>(prefix_popcount(k));
    const long double denom = static_cast<long double>(1ULL << level);
    return numer / denom;
}

long double overlap_height(const long double x, const long double curve_y) {
    const long double dx = x - 0.25L;
    const long double inside = 0.0625L - dx * dx;
    if (inside <= 0.0L) {
        return 0.0L;
    }

    const long double dy = std::sqrt(inside);
    const long double circle_low = 0.5L - dy;
    const long double circle_high = 0.5L + dy;
    const long double lo = std::max(0.0L, circle_low);
    const long double hi = std::min(circle_high, curve_y);
    return hi > lo ? hi - lo : 0.0L;
}

struct PartialSums {
    long double odd_sum = 0.0L;
    long double even_sum = 0.0L;
};

PartialSums integrate_chunk(const int level, const u64 begin, const u64 end) {
    PartialSums out;
    if (begin >= end) {
        return out;
    }

    const long double inv_denom = 1.0L / static_cast<long double>(1ULL << level);
    long double numer = static_cast<long double>(level) * static_cast<long double>(begin) -
                        2.0L * static_cast<long double>(prefix_popcount(begin));

    for (u64 i = begin; i < end; ++i) {
        const long double x = static_cast<long double>(i) * inv_denom;
        const long double curve_y = numer * inv_denom;
        const long double fx = overlap_height(x, curve_y);
        if ((i & 1ULL) != 0ULL) {
            out.odd_sum += fx;
        } else {
            out.even_sum += fx;
        }
        numer += static_cast<long double>(level - 2 * static_cast<int>(__builtin_popcountll(i)));
    }

    return out;
}

long double simpson_level(const int level, int threads) {
    const u64 denom = 1ULL << level;
    const u64 steps = denom >> 1U;
    if (steps < 2) {
        return 0.0L;
    }

    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
        if (threads <= 0) {
            threads = 1;
        }
    }

    const u64 interior = steps - 1;
    if (interior == 0) {
        return 0.0L;
    }
    threads = std::max(1, std::min<int>(threads, static_cast<int>(interior)));

    std::vector<PartialSums> partials(static_cast<std::size_t>(threads));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));

    const u64 chunk = (interior + static_cast<u64>(threads) - 1ULL) / static_cast<u64>(threads);
    for (int t = 0; t < threads; ++t) {
        const u64 begin = 1ULL + static_cast<u64>(t) * chunk;
        const u64 end = std::min(steps, begin + chunk);
        if (begin >= end) {
            continue;
        }
        workers.emplace_back([&, t, begin, end]() { partials[static_cast<std::size_t>(t)] = integrate_chunk(level, begin, end); });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    long double odd_sum = 0.0L;
    long double even_sum = 0.0L;
    for (const PartialSums& p : partials) {
        odd_sum += p.odd_sum;
        even_sum += p.even_sum;
    }

    const long double h = 1.0L / static_cast<long double>(denom);
    const long double f0 = 0.0L;
    const long double fn = overlap_height(0.5L, 0.5L);
    return h * (f0 + fn + 4.0L * odd_sum + 2.0L * even_sum) / 3.0L;
}

long double solve(const long double eps, const int min_level, const int max_level, const int threads) {
    long double previous = simpson_level(min_level, threads);
    for (int level = min_level + 1; level <= max_level; ++level) {
        const long double current = simpson_level(level, threads);
        if (std::fabsl(current - previous) < eps) {
            return current;
        }
        previous = current;
    }
    return previous;
}

bool run_checkpoints() {
    if (prefix_popcount(8) != 12 || prefix_popcount(16) != 32) {
        std::cerr << "Checkpoint failed for prefix popcount" << '\n';
        return false;
    }
    if (std::fabsl(takagi_dyadic(1, 1) - 0.5L) > 1e-15L || std::fabsl(takagi_dyadic(1, 2) - 0.5L) > 1e-15L ||
        std::fabsl(takagi_dyadic(1, 3) - 0.375L) > 1e-15L || std::fabsl(takagi_dyadic(3, 3) - 0.625L) > 1e-15L) {
        std::cerr << "Checkpoint failed for dyadic Takagi values" << '\n';
        return false;
    }

    const long double coarse = simpson_level(20, 1);
    const long double fine = simpson_level(21, 1);
    if (coarse <= 0.0L || fine <= 0.0L || fine >= 0.2L || std::fabsl(fine - coarse) > 1e-4L) {
        std::cerr << "Checkpoint failed for area refinement" << '\n';
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

    const long double answer = solve(options.eps, options.min_level, options.max_level, options.threads);
    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(answer) << '\n';
    return 0;
}
