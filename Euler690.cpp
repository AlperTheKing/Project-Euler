#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;

constexpr i64 kMod = 1'000'000'007LL;
constexpr int kDefaultN = 2019;

struct Options {
    int n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

inline i64 add_mod(const i64 a, const i64 b) {
    i64 out = a + b;
    if (out >= kMod) {
        out -= kMod;
    }
    return out;
}

inline i64 sub_mod(const i64 a, const i64 b) {
    i64 out = a - b;
    if (out < 0) {
        out += kMod;
    }
    return out;
}

inline i64 mul_mod(const i64 a, const i64 b) {
    return static_cast<i64>((static_cast<__int128>(a) * static_cast<__int128>(b)) % kMod);
}

i64 pow_mod(i64 base, i64 exp) {
    i64 out = 1;
    while (exp > 0) {
        if ((exp & 1LL) != 0LL) {
            out = mul_mod(out, base);
        }
        base = mul_mod(base, base);
        exp >>= 1LL;
    }
    return out;
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::uint64_t parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
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

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--n=", parsed_u32)) {
            if (parsed_u32 > static_cast<u32>(std::numeric_limits<int>::max())) {
                std::cerr << "--n is too large.\n";
                return false;
            }
            options.n = static_cast<int>(parsed_u32);
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

    if (options.n < 1) {
        std::cerr << "--n must be at least 1.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const int workload) {
    if (!allow_multithreading || workload < 512) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

std::vector<i64> precompute_inverses(const int n) {
    std::vector<i64> inv(static_cast<std::size_t>(n) + 1ULL, 0LL);
    if (n >= 1) {
        inv[1] = 1LL;
    }

    for (int i = 2; i <= n; ++i) {
        const i64 quotient = kMod / static_cast<i64>(i);
        const i64 remainder = kMod % static_cast<i64>(i);
        inv[static_cast<std::size_t>(i)] =
            sub_mod(0LL, mul_mod(quotient, inv[static_cast<std::size_t>(remainder)]));
    }

    return inv;
}

std::vector<i64> convolve_truncated(const std::vector<i64>& a,
                                    const std::vector<i64>& b,
                                    const int n,
                                    const bool allow_multithreading,
                                    const unsigned requested_threads) {
    std::vector<i64> out(static_cast<std::size_t>(n) + 1ULL, 0LL);
    if (a.empty() || b.empty() || n < 0) {
        return out;
    }

    const int max_i = std::min<int>(n, static_cast<int>(a.size()) - 1);
    const int max_j_total = std::min<int>(n, static_cast<int>(b.size()) - 1);
    if (max_i < 0 || max_j_total < 0) {
        return out;
    }

    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, max_i + 1);
    if (threads == 1U) {
        for (int i = 0; i <= max_i; ++i) {
            const i64 ai = a[static_cast<std::size_t>(i)];
            if (ai == 0LL) {
                continue;
            }
            const int max_j = std::min(max_j_total, n - i);
            for (int j = 0; j <= max_j; ++j) {
                const i64 bj = b[static_cast<std::size_t>(j)];
                if (bj == 0LL) {
                    continue;
                }
                const std::size_t idx = static_cast<std::size_t>(i + j);
                out[idx] = add_mod(out[idx], mul_mod(ai, bj));
            }
        }
        return out;
    }

    std::vector<std::vector<i64>> partial(
        threads, std::vector<i64>(static_cast<std::size_t>(n) + 1ULL, 0LL));
    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        const int start = static_cast<int>((static_cast<std::int64_t>(t) * (max_i + 1)) / threads);
        const int end =
            static_cast<int>((static_cast<std::int64_t>(t + 1U) * (max_i + 1)) / threads);

        pool.emplace_back([&, t, start, end]() {
            std::vector<i64>& local = partial[static_cast<std::size_t>(t)];
            for (int i = start; i < end; ++i) {
                const i64 ai = a[static_cast<std::size_t>(i)];
                if (ai == 0LL) {
                    continue;
                }
                const int max_j = std::min(max_j_total, n - i);
                for (int j = 0; j <= max_j; ++j) {
                    const i64 bj = b[static_cast<std::size_t>(j)];
                    if (bj == 0LL) {
                        continue;
                    }
                    const std::size_t idx = static_cast<std::size_t>(i + j);
                    local[idx] = add_mod(local[idx], mul_mod(ai, bj));
                }
            }
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    for (unsigned t = 0; t < threads; ++t) {
        const std::vector<i64>& local = partial[static_cast<std::size_t>(t)];
        for (int i = 0; i <= n; ++i) {
            out[static_cast<std::size_t>(i)] =
                add_mod(out[static_cast<std::size_t>(i)], local[static_cast<std::size_t>(i)]);
        }
    }

    return out;
}

std::vector<i64> invert_series(const std::vector<i64>& den, const int n) {
    std::vector<i64> out(static_cast<std::size_t>(n) + 1ULL, 0LL);
    out[0] = pow_mod(den[0], kMod - 2);

    for (int i = 1; i <= n; ++i) {
        i64 sum = 0LL;
        for (int k = 1; k <= i; ++k) {
            sum = add_mod(sum, mul_mod(den[static_cast<std::size_t>(k)],
                                       out[static_cast<std::size_t>(i - k)]));
        }
        out[static_cast<std::size_t>(i)] = sub_mod(0LL, mul_mod(out[0], sum));
    }

    return out;
}

std::vector<i64> partition_series(const int n) {
    std::vector<i64> p(static_cast<std::size_t>(n) + 1ULL, 0LL);
    p[0] = 1LL;

    for (int part = 1; part <= n; ++part) {
        for (int i = part; i <= n; ++i) {
            p[static_cast<std::size_t>(i)] =
                add_mod(p[static_cast<std::size_t>(i)], p[static_cast<std::size_t>(i - part)]);
        }
    }

    return p;
}

std::vector<i64> compute_lobster_tree_counts(const int n,
                                             const std::vector<i64>& inv,
                                             const bool allow_multithreading,
                                             const unsigned requested_threads) {
    const std::vector<i64> p = partition_series(n);

    std::vector<i64> inv_1mx(static_cast<std::size_t>(n) + 1ULL, 1LL);
    std::vector<i64> a(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int i = 0; i <= n; ++i) {
        a[static_cast<std::size_t>(i)] =
            sub_mod(p[static_cast<std::size_t>(i)], inv_1mx[static_cast<std::size_t>(i)]);
    }

    std::vector<i64> x_p(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int i = 1; i <= n; ++i) {
        x_p[static_cast<std::size_t>(i)] = p[static_cast<std::size_t>(i - 1)];
    }

    std::vector<i64> den1(static_cast<std::size_t>(n) + 1ULL, 0LL);
    den1[0] = 1LL;
    for (int i = 1; i <= n; ++i) {
        den1[static_cast<std::size_t>(i)] = sub_mod(0LL, x_p[static_cast<std::size_t>(i)]);
    }
    const std::vector<i64> inv_den1 = invert_series(den1, n);

    const std::vector<i64> a_sq =
        convolve_truncated(a, a, n, allow_multithreading, requested_threads);
    const std::vector<i64> part1 =
        convolve_truncated(a_sq, inv_den1, n, allow_multithreading, requested_threads);

    std::vector<i64> p2(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int i = 0; i <= n; i += 2) {
        p2[static_cast<std::size_t>(i)] = p[static_cast<std::size_t>(i / 2)];
    }

    std::vector<i64> inv_1mx2(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int i = 0; i <= n; i += 2) {
        inv_1mx2[static_cast<std::size_t>(i)] = 1LL;
    }

    std::vector<i64> b(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int i = 0; i <= n; ++i) {
        b[static_cast<std::size_t>(i)] =
            sub_mod(p2[static_cast<std::size_t>(i)], inv_1mx2[static_cast<std::size_t>(i)]);
    }

    std::vector<i64> one_plus_x_p = x_p;
    one_plus_x_p[0] = add_mod(one_plus_x_p[0], 1LL);

    std::vector<i64> den2(static_cast<std::size_t>(n) + 1ULL, 0LL);
    den2[0] = 1LL;
    for (int i = 2; i <= n; ++i) {
        den2[static_cast<std::size_t>(i)] = sub_mod(0LL, p2[static_cast<std::size_t>(i - 2)]);
    }
    const std::vector<i64> inv_den2 = invert_series(den2, n);

    const std::vector<i64> tmp =
        convolve_truncated(b, one_plus_x_p, n, allow_multithreading, requested_threads);
    const std::vector<i64> part2 =
        convolve_truncated(tmp, inv_den2, n, allow_multithreading, requested_threads);

    std::vector<i64> lobsters(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int i = 0; i + 2 <= n; ++i) {
        const i64 sum = add_mod(part1[static_cast<std::size_t>(i)], part2[static_cast<std::size_t>(i)]);
        lobsters[static_cast<std::size_t>(i + 2)] =
            add_mod(lobsters[static_cast<std::size_t>(i + 2)],
                    mul_mod(sum, inv[2]));
    }

    for (int i = 1; i <= n; ++i) {
        lobsters[static_cast<std::size_t>(i)] =
            add_mod(lobsters[static_cast<std::size_t>(i)], x_p[static_cast<std::size_t>(i)]);
    }

    std::vector<i64> rec(static_cast<std::size_t>(n) + 1ULL, 0LL);
    rec[0] = 1LL;
    for (int i = 1; i <= n; ++i) {
        i64 v = rec[static_cast<std::size_t>(i - 1)];
        if (i >= 2) {
            v = add_mod(v, rec[static_cast<std::size_t>(i - 2)]);
        }
        if (i >= 3) {
            v = sub_mod(v, rec[static_cast<std::size_t>(i - 3)]);
        }
        rec[static_cast<std::size_t>(i)] = v;
    }

    for (int i = 0; i + 3 <= n; ++i) {
        lobsters[static_cast<std::size_t>(i + 3)] =
            sub_mod(lobsters[static_cast<std::size_t>(i + 3)], rec[static_cast<std::size_t>(i)]);
    }

    return lobsters;
}

std::vector<i64> euler_transform_from_tree_counts(const std::vector<i64>& tree_counts,
                                                   const int n,
                                                   const std::vector<i64>& inv) {
    std::vector<i64> c(static_cast<std::size_t>(n) + 1ULL, 0LL);
    for (int d = 1; d <= n; ++d) {
        const i64 val = mul_mod(static_cast<i64>(d), tree_counts[static_cast<std::size_t>(d)]);
        for (int m = d; m <= n; m += d) {
            c[static_cast<std::size_t>(m)] = add_mod(c[static_cast<std::size_t>(m)], val);
        }
    }

    std::vector<i64> out(static_cast<std::size_t>(n) + 1ULL, 0LL);
    out[0] = 1LL;

    for (int i = 1; i <= n; ++i) {
        i64 sum = 0LL;
        for (int k = 1; k <= i; ++k) {
            sum = add_mod(sum,
                          mul_mod(c[static_cast<std::size_t>(k)],
                                  out[static_cast<std::size_t>(i - k)]));
        }
        out[static_cast<std::size_t>(i)] = mul_mod(sum, inv[static_cast<std::size_t>(i)]);
    }

    return out;
}

bool run_checkpoints() {
    constexpr int kCheckN = 20;
    const std::vector<i64> inv = precompute_inverses(kCheckN);
    const std::vector<i64> lobster_trees =
        compute_lobster_tree_counts(kCheckN, inv, false, 1U);

    const std::vector<i64> expected_tree_prefix = {
        0LL, 1LL, 1LL, 1LL, 2LL, 3LL, 6LL, 11LL, 23LL, 47LL, 105LL, 231LL, 532LL};
    for (std::size_t i = 1; i < expected_tree_prefix.size(); ++i) {
        if (lobster_trees[i] != expected_tree_prefix[i]) {
            std::cerr << "Checkpoint failed: lobster tree count mismatch at n=" << i << ".\n";
            return false;
        }
    }

    const std::vector<i64> tom = euler_transform_from_tree_counts(lobster_trees, kCheckN, inv);

    if (tom[3] != 3LL) {
        std::cerr << "Checkpoint failed: T(3) != 3.\n";
        return false;
    }
    if (tom[7] != 37LL) {
        std::cerr << "Checkpoint failed: T(7) != 37.\n";
        return false;
    }
    if (tom[10] != 328LL) {
        std::cerr << "Checkpoint failed: T(10) != 328.\n";
        return false;
    }
    if (tom[20] != 1'416'269LL) {
        std::cerr << "Checkpoint failed: T(20) != 1416269.\n";
        return false;
    }

    return true;
}

i64 solve(const int n, const bool allow_multithreading, const unsigned requested_threads) {
    const std::vector<i64> inv = precompute_inverses(n);
    const std::vector<i64> lobster_trees =
        compute_lobster_tree_counts(n, inv, allow_multithreading, requested_threads);
    const std::vector<i64> tom = euler_transform_from_tree_counts(lobster_trees, n, inv);
    return tom[static_cast<std::size_t>(n)];
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }

    const i64 answer = solve(options.n, options.allow_multithreading, options.requested_threads);
    std::cout << answer << '\n';
    return 0;
}
