#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

struct Options {
    i64 limit = 100000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 3;
}

u64 gcd3(const u64 a, const u64 b, const u64 c) {
    return std::gcd(a, std::gcd(b, c));
}

struct Counts {
    u64 with_60 = 0;
    u64 with_90 = 0;
    u64 with_120 = 0;

    u64 total() const {
        return with_60 + with_90 + with_120;
    }
};

u64 count_90_range(const i64 limit, const i64 m_start, const i64 m_end, const i64 m_step) {
    u64 out = 0;
    for (i64 m = m_start; m <= m_end; m += m_step) {
        for (i64 n = 1; n < m; ++n) {
            if (((m - n) & 1LL) == 0 || std::gcd(m, n) != 1) {
                continue;
            }
            const i64 primitive_perimeter = 2 * m * (m + n);
            out += static_cast<u64>(limit / primitive_perimeter);
        }
    }
    return out;
}

u64 count_120_range(const i64 limit, const i64 m_start, const i64 m_end, const i64 m_step) {
    u64 out = 0;
    for (i64 m = m_start; m <= m_end; m += m_step) {
        for (i64 n = 1; n < m; ++n) {
            if (std::gcd(m, n) != 1) {
                continue;
            }
            const i64 mm = m * m;
            const i64 nn = n * n;
            const i64 a = mm - nn;
            const i64 b = 2 * m * n + nn;
            const i64 c = mm + m * n + nn;
            if (gcd3(static_cast<u64>(a), static_cast<u64>(b), static_cast<u64>(c)) != 1) {
                continue;
            }
            const i64 primitive_perimeter = a + b + c;
            out += static_cast<u64>(limit / primitive_perimeter);
        }
    }
    return out;
}

u64 count_60_non_eq_range(const i64 limit, const i64 m_start, const i64 m_end, const i64 m_step) {
    u64 out = 0;
    for (i64 m = m_start; m <= m_end; m += m_step) {
        const i64 mm = m * m;
        const i64 n_max = (m - 1) / 2;
        for (i64 n = 1; n <= n_max; ++n) {
            if (std::gcd(m, n) != 1) {
                continue;
            }
            const i64 nn = n * n;
            const i64 a = mm - nn;
            const i64 c = mm - m * n + nn;

            const i64 b_a = 2 * m * n - nn;
            if (gcd3(static_cast<u64>(a), static_cast<u64>(b_a), static_cast<u64>(c)) == 1) {
                const i64 primitive_perimeter = a + b_a + c;
                out += static_cast<u64>(limit / primitive_perimeter);
            }

            const i64 b_b = mm - 2 * m * n;
            if (b_b > 0 && gcd3(static_cast<u64>(a), static_cast<u64>(b_b), static_cast<u64>(c)) == 1) {
                const i64 primitive_perimeter = a + b_b + c;
                out += static_cast<u64>(limit / primitive_perimeter);
            }
        }
    }
    return out;
}

i64 max_m_90(const i64 limit) {
    i64 m = 2;
    while (2 * m * (m + 1) <= limit) {
        ++m;
    }
    return m - 1;
}

i64 max_m_120(const i64 limit) {
    i64 m = 2;
    while (2 * m * m + 3 * m + 1 <= limit) {
        ++m;
    }
    return m - 1;
}

i64 max_m_60(const i64 limit) {
    i64 m = 2;
    while (true) {
        const __int128 mm = static_cast<__int128>(m) * static_cast<__int128>(m);
        const __int128 min_perimeter_a = 2 * mm + m - 1;
        const __int128 min_perimeter_b = static_cast<__int128>(3) * m * ((m + 1) / 2);
        if (min_perimeter_a > limit && min_perimeter_b > limit) {
            break;
        }
        ++m;
    }
    return m - 1;
}

struct WorkerCtx {
    i64 limit = 0;
    i64 m_max_90 = 0;
    i64 m_max_120 = 0;
    i64 m_max_60 = 0;
    i64 m_start = 0;
    i64 m_step = 0;
    u64 c90 = 0;
    u64 c120 = 0;
    u64 c60 = 0;
};

void* worker_main(void* ptr) {
    auto* ctx = static_cast<WorkerCtx*>(ptr);
    ctx->c90 = count_90_range(ctx->limit, ctx->m_start, ctx->m_max_90, ctx->m_step);
    ctx->c120 = count_120_range(ctx->limit, ctx->m_start, ctx->m_max_120, ctx->m_step);
    ctx->c60 = count_60_non_eq_range(ctx->limit, ctx->m_start, ctx->m_max_60, ctx->m_step);
    return nullptr;
}

unsigned choose_thread_count() {
    long cpu = sysconf(_SC_NPROCESSORS_ONLN);
    unsigned threads = (cpu > 0) ? static_cast<unsigned>(cpu) : 1U;
    if (threads > 8U) {
        threads = 8U;
    }
    if (threads == 0U) {
        threads = 1U;
    }
    return threads;
}

Counts count_by_parameterization(const i64 limit) {
    Counts counts;
    const i64 m_max_90 = max_m_90(limit);
    const i64 m_max_120 = max_m_120(limit);
    const i64 m_max_60 = max_m_60(limit);

    counts.with_60 = static_cast<u64>(limit / 3);
    const unsigned thread_count = choose_thread_count();
    if (thread_count <= 1U) {
        counts.with_90 = count_90_range(limit, 2, m_max_90, 1);
        counts.with_120 = count_120_range(limit, 2, m_max_120, 1);
        counts.with_60 += count_60_non_eq_range(limit, 2, m_max_60, 1);
        return counts;
    }

    std::vector<pthread_t> threads(thread_count);
    std::vector<WorkerCtx> ctx(thread_count);
    unsigned created = 0U;
    bool failed = false;
    for (unsigned t = 0; t < thread_count; ++t) {
        ctx[t] = WorkerCtx{limit,
                           m_max_90,
                           m_max_120,
                           m_max_60,
                           static_cast<i64>(2 + t),
                           static_cast<i64>(thread_count),
                           0ULL,
                           0ULL,
                           0ULL};
        if (pthread_create(&threads[t], nullptr, worker_main, &ctx[t]) != 0) {
            failed = true;
            break;
        }
        ++created;
    }
    for (unsigned t = 0; t < created; ++t) {
        pthread_join(threads[t], nullptr);
    }

    if (failed) {
        counts.with_90 = count_90_range(limit, 2, m_max_90, 1);
        counts.with_120 = count_120_range(limit, 2, m_max_120, 1);
        counts.with_60 += count_60_non_eq_range(limit, 2, m_max_60, 1);
        return counts;
    }

    for (const WorkerCtx& w : ctx) {
        counts.with_90 += w.c90;
        counts.with_120 += w.c120;
        counts.with_60 += w.c60;
    }

    return counts;
}

Counts count_by_bruteforce(const int limit) {
    Counts counts;

    for (int a = 1; a <= limit / 3; ++a) {
        for (int b = a; b <= (limit - a) / 2; ++b) {
            const int c_max = limit - a - b;
            for (int c = b; c <= c_max; ++c) {
                if (a + b <= c) {
                    continue;
                }

                const std::array<i64, 3> sides = {a, b, c};
                bool has_60 = false;
                bool has_90 = false;
                bool has_120 = false;

                for (int i = 0; i < 3; ++i) {
                    const i64 z = sides[static_cast<std::size_t>(i)];
                    const i64 x = sides[static_cast<std::size_t>((i + 1) % 3)];
                    const i64 y = sides[static_cast<std::size_t>((i + 2) % 3)];
                    const i64 zz = z * z;
                    const i64 xx_yy = x * x + y * y;

                    if (zz == xx_yy) {
                        has_90 = true;
                    }
                    if (zz == xx_yy - x * y) {
                        has_60 = true;
                    }
                    if (zz == xx_yy + x * y) {
                        has_120 = true;
                    }
                }

                counts.with_60 += static_cast<u64>(has_60);
                counts.with_90 += static_cast<u64>(has_90);
                counts.with_120 += static_cast<u64>(has_120);
            }
        }
    }

    return counts;
}

bool run_checkpoints() {
    {
        const Counts fast = count_by_parameterization(100);
        if (fast.total() != 84 || fast.with_60 != 53 || fast.with_90 != 17 || fast.with_120 != 14) {
            std::cerr << "Checkpoint failed at limit=100" << '\n';
            return false;
        }
    }

    for (const int limit : {120, 200, 300}) {
        const Counts brute = count_by_bruteforce(limit);
        const Counts fast = count_by_parameterization(limit);
        if (brute.with_60 != fast.with_60 || brute.with_90 != fast.with_90 ||
            brute.with_120 != fast.with_120 || brute.total() != fast.total()) {
            std::cerr << "Mismatch at limit=" << limit << ": brute(total=" << brute.total()
                      << ", 60=" << brute.with_60 << ", 90=" << brute.with_90
                      << ", 120=" << brute.with_120 << ") fast(total=" << fast.total()
                      << ", 60=" << fast.with_60 << ", 90=" << fast.with_90
                      << ", 120=" << fast.with_120 << ")" << '\n';
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

    const Counts result = count_by_parameterization(options.limit);
    std::cout << result.total() << '\n';
    return 0;
}
