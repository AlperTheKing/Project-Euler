#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kDefaultN = 100'000'000ULL;
constexpr u64 kDefaultDirectLimit = 220'000ULL;

struct Options {
    u64 n = kDefaultN;
    u64 direct_limit = kDefaultDirectLimit;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
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
            options.n = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--direct-limit=", parsed_u64)) {
            options.direct_limit = parsed_u64;
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

    if (options.direct_limit < 2ULL) {
        std::cerr << "--direct-limit must be >= 2.\n";
        return false;
    }
    if (options.direct_limit > static_cast<u64>(std::numeric_limits<int>::max())) {
        std::cerr << "--direct-limit is too large for this implementation.\n";
        return false;
    }
    if (options.run_checkpoints && options.direct_limit < 10'000ULL) {
        std::cerr << "--direct-limit must be >= 10000 when checkpoints are enabled.\n";
        return false;
    }

    return true;
}

u64 next_sequence_value(const u64 x) {
    return static_cast<u64>((static_cast<u128>(x) * static_cast<u128>(x) + 45ULL) % kMod);
}

struct CycleInfo {
    u64 mu = 0ULL;
    u64 lambda = 0ULL;
};

CycleInfo find_cycle_brent(const u64 x0) {
    u64 power = 1ULL;
    u64 lambda = 1ULL;
    u64 tortoise = x0;
    u64 hare = next_sequence_value(x0);

    while (tortoise != hare) {
        if (power == lambda) {
            tortoise = hare;
            power <<= 1U;
            lambda = 0ULL;
        }
        hare = next_sequence_value(hare);
        ++lambda;
    }

    u64 mu = 0ULL;
    tortoise = x0;
    hare = x0;
    for (u64 i = 0ULL; i < lambda; ++i) {
        hare = next_sequence_value(hare);
    }
    while (tortoise != hare) {
        tortoise = next_sequence_value(tortoise);
        hare = next_sequence_value(hare);
        ++mu;
    }

    return CycleInfo{mu, lambda};
}

std::vector<std::uint32_t> generate_sequence_prefix(const int n) {
    std::vector<std::uint32_t> seq(static_cast<std::size_t>(n), 0U);
    if (n <= 1) {
        return seq;
    }

    u64 x = 0ULL;
    for (int i = 1; i < n; ++i) {
        x = next_sequence_value(x);
        seq[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(x);
    }
    return seq;
}

std::vector<i64> compute_f_table(const std::vector<std::uint32_t>& seq) {
    const int n = static_cast<int>(seq.size());
    std::vector<i64> f(static_cast<std::size_t>(n + 1), 0LL);
    if (n == 0) {
        return f;
    }

    std::vector<i64> prefix(static_cast<std::size_t>(n + 1), 0LL);
    std::vector<i64> prev(static_cast<std::size_t>(n), 0LL);
    std::vector<i64> cur(static_cast<std::size_t>(n), 0LL);

    for (int i = 0; i < n; ++i) {
        prefix[static_cast<std::size_t>(i + 1)] =
            prefix[static_cast<std::size_t>(i)] + static_cast<i64>(seq[static_cast<std::size_t>(i)]);
        prev[static_cast<std::size_t>(i)] = static_cast<i64>(seq[static_cast<std::size_t>(i)]);
    }

    f[1] = (prefix[1] + prev[0]) / 2LL;

    for (int len = 2; len <= n; ++len) {
        const int upto = n - len;
        int right = len - 1;

        for (int i = 0; i <= upto; ++i, ++right) {
            const i64 left_pick =
                static_cast<i64>(seq[static_cast<std::size_t>(i)]) -
                prev[static_cast<std::size_t>(i + 1)];
            const i64 right_pick =
                static_cast<i64>(seq[static_cast<std::size_t>(right)]) -
                prev[static_cast<std::size_t>(i)];
            cur[static_cast<std::size_t>(i)] = (left_pick > right_pick) ? left_pick : right_pick;
        }

        std::swap(prev, cur);
        f[static_cast<std::size_t>(len)] =
            (prefix[static_cast<std::size_t>(len)] + prev[0]) / 2LL;
    }

    return f;
}

bool run_problem_checkpoints(const std::vector<i64>& f) {
    struct Checkpoint {
        int n;
        i64 expected;
    };

    const Checkpoint checkpoints[] = {
        {2, 45LL},
        {4, 4'284'990LL},
        {100, 26'365'463'243LL},
        {10'000, 2'495'838'522'951LL},
    };

    for (const auto& checkpoint : checkpoints) {
        if (checkpoint.n >= static_cast<int>(f.size())) {
            std::cerr << "Checkpoint n=" << checkpoint.n
                      << " is out of range for the current --direct-limit.\n";
            return false;
        }
        const i64 got = f[static_cast<std::size_t>(checkpoint.n)];
        if (got != checkpoint.expected) {
            std::cerr << "Checkpoint failed at n=" << checkpoint.n
                      << ": got " << got << ", expected " << checkpoint.expected << '\n';
            return false;
        }
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading, const unsigned requested_threads) {
    if (!allow_multithreading) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, threads);
}

struct AffineInfo {
    u64 start = 0ULL;
    i64 delta = 0LL;
    bool valid = false;
};

AffineInfo find_affine_for_parity(const std::vector<i64>& f,
                                  const u64 period,
                                  const int parity) {
    const u64 nmax = static_cast<u64>(f.size() - 1U);
    const u64 first = (parity == 0) ? 2ULL : 1ULL;  // parity=0 -> even n, parity=1 -> odd n

    if (first + period > nmax) {
        return AffineInfo{};
    }

    for (u64 start = first; start + period <= nmax; start += 2ULL) {
        const i64 delta =
            f[static_cast<std::size_t>(start + period)] - f[static_cast<std::size_t>(start)];
        bool ok = true;
        for (u64 n = start; n + period <= nmax; n += 2ULL) {
            const i64 got =
                f[static_cast<std::size_t>(n + period)] - f[static_cast<std::size_t>(n)];
            if (got != delta) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return AffineInfo{start, delta, true};
        }
    }

    return AffineInfo{};
}

std::array<AffineInfo, 2> find_affine_infos(const std::vector<i64>& f,
                                            const u64 period,
                                            const bool allow_multithreading,
                                            const unsigned requested_threads) {
    std::array<AffineInfo, 2> infos{};
    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads);

    if (threads >= 2U) {
        std::thread even_worker([&]() {
            infos[0] = find_affine_for_parity(f, period, 0);  // even n
        });
        infos[1] = find_affine_for_parity(f, period, 1);      // odd n
        even_worker.join();
    } else {
        infos[0] = find_affine_for_parity(f, period, 0);
        infos[1] = find_affine_for_parity(f, period, 1);
    }

    return infos;
}

i64 extrapolated_f(const u64 n,
                   const std::vector<i64>& f,
                   const u64 period,
                   const std::array<AffineInfo, 2>& infos) {
    if (n < f.size()) {
        return f[static_cast<std::size_t>(n)];
    }

    const int parity = static_cast<int>(n & 1ULL);  // 0 even, 1 odd
    const AffineInfo& info = infos[parity];
    if (!info.valid || n < info.start) {
        return -1LL;
    }

    i64 rem = static_cast<i64>((n - info.start) % period);
    if (rem < 0LL) {
        rem += static_cast<i64>(period);
    }

    const u64 base = info.start + static_cast<u64>(rem);
    const u64 blocks = (n - base) / period;
    return f[static_cast<std::size_t>(base)] + static_cast<i64>(blocks) * info.delta;
}

bool validate_affine_model(const std::vector<i64>& f,
                           const u64 period,
                           const std::array<AffineInfo, 2>& infos) {
    const u64 nmax = static_cast<u64>(f.size() - 1U);
    for (int parity = 0; parity < 2; ++parity) {
        if (!infos[parity].valid) {
            std::cerr << "Could not detect affine-periodic behavior for parity " << parity
                      << ". Increase --direct-limit.\n";
            return false;
        }
        for (u64 n = infos[parity].start; n + period <= nmax; n += 2ULL) {
            const i64 got =
                f[static_cast<std::size_t>(n + period)] - f[static_cast<std::size_t>(n)];
            if (got != infos[parity].delta) {
                std::cerr << "Affine validation failed at n=" << n
                          << " for parity " << parity << ".\n";
                return false;
            }
        }
    }

    // Extra tail check: the extrapolator must match direct DP values near the window end.
    const u64 tail_window = std::min<u64>(nmax, 5ULL * period);
    const u64 begin = nmax - tail_window;
    for (u64 n = begin; n <= nmax; ++n) {
        const i64 predicted = extrapolated_f(n, f, period, infos);
        if (predicted != f[static_cast<std::size_t>(n)]) {
            std::cerr << "Tail extrapolation mismatch at n=" << n
                      << ": got " << predicted
                      << ", expected " << f[static_cast<std::size_t>(n)] << ".\n";
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

    const CycleInfo cycle = find_cycle_brent(0ULL);
    if (cycle.lambda == 0ULL) {
        std::cerr << "Cycle detection failed.\n";
        return 1;
    }

    const u64 period = 2ULL * cycle.lambda;
    if (period > options.direct_limit) {
        std::cerr << "--direct-limit must be at least " << period
                  << " for period detection.\n";
        return 1;
    }

    const int direct_n = static_cast<int>(options.direct_limit);
    const std::vector<std::uint32_t> seq = generate_sequence_prefix(direct_n);
    const std::vector<i64> f = compute_f_table(seq);

    if (options.run_checkpoints) {
        if (cycle.mu != 57'956ULL || cycle.lambda != 7'248ULL) {
            std::cerr << "Cycle checkpoint failed: mu=" << cycle.mu
                      << ", lambda=" << cycle.lambda << ".\n";
            return 1;
        }
        if (!run_problem_checkpoints(f)) {
            return 1;
        }
    }

    const auto infos =
        find_affine_infos(f, period, options.allow_multithreading, options.requested_threads);
    if (!validate_affine_model(f, period, infos)) {
        return 1;
    }

    const i64 answer = extrapolated_f(options.n, f, period, infos);
    if (answer < 0LL) {
        std::cerr << "Could not extrapolate F(" << options.n
                  << "). Try increasing --direct-limit.\n";
        return 1;
    }

    std::cout << answer << '\n';
    return 0;
}
