#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;

constexpr unsigned kDefaultN = 14U;
constexpr unsigned kMaxPracticalN = 16U;

struct Options {
    unsigned n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct CompetitionResult {
    std::vector<long double> win_from_start;
    long double expected_dishes_from_start = 0.0L;
};

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u32 parsed = 0U;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u32 digit = static_cast<u32>(ch - '0');
        if (parsed > (std::numeric_limits<u32>::max() - digit) / 10U) {
            return false;
        }
        parsed = parsed * 10U + digit;
    }

    value = parsed;
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

        u32 parsed = 0U;
        if (parse_u32_after_prefix(arg, "--n=", parsed)) {
            options.n = parsed;
            continue;
        }
        if (parse_u32_after_prefix(arg, "--threads=", parsed)) {
            options.requested_threads = parsed;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n < 2U) {
        std::cerr << "--n must be >= 2.\n";
        return false;
    }
    if (options.n > kMaxPracticalN) {
        std::cerr << "--n too large for this exact implementation (max "
                  << kMaxPracticalN << ").\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t jobs) {
    constexpr std::size_t kMinJobsForParallel = 128U;
    if (!allow_multithreading || jobs < kMinJobsForParallel) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(jobs)));
}

std::vector<long double> fibonacci_skills(const unsigned n) {
    std::vector<long double> fib(static_cast<std::size_t>(n + 2U), 0.0L);
    fib[1] = 1.0L;
    fib[2] = 1.0L;
    for (unsigned i = 3U; i <= n + 1U; ++i) {
        fib[i] = fib[i - 1U] + fib[i - 2U];
    }

    std::vector<long double> skills(static_cast<std::size_t>(n), 0.0L);
    for (unsigned k = 1U; k <= n; ++k) {
        skills[static_cast<std::size_t>(k - 1U)] =
            fib[k] / fib[n + 1U];
    }
    return skills;
}

inline std::size_t state_index(const u32 mask, const unsigned turn, const unsigned n) {
    return static_cast<std::size_t>(mask) * n + turn;
}

inline unsigned next_alive_after(const u32 mask, const unsigned chef) {
    const u32 higher_bits = mask & (~((1U << (chef + 1U)) - 1U));
    if (higher_bits != 0U) {
        return static_cast<unsigned>(__builtin_ctz(higher_bits));
    }
    return static_cast<unsigned>(__builtin_ctz(mask));
}

CompetitionResult solve_competition(const std::vector<long double>& skills,
                                    const bool allow_multithreading,
                                    const unsigned requested_threads) {
    const unsigned n = static_cast<unsigned>(skills.size());
    const u32 total_masks = (1U << n);
    const std::size_t total_states = static_cast<std::size_t>(total_masks) * n;
    const std::size_t win_stride = n;
    const std::size_t win_size = total_states * win_stride;

    std::vector<long double> wins(win_size, 0.0L);
    std::vector<long double> dishes(total_states, 0.0L);
    std::vector<std::vector<u32>> masks_by_size(static_cast<std::size_t>(n + 1U));

    for (u32 mask = 1U; mask < total_masks; ++mask) {
        const unsigned bits = static_cast<unsigned>(__builtin_popcount(mask));
        masks_by_size[bits].push_back(mask);
    }

    const auto process_mask = [&](const u32 mask) {
        std::vector<unsigned> members;
        members.reserve(n);

        u32 temp = mask;
        while (temp != 0U) {
            const unsigned chef = static_cast<unsigned>(__builtin_ctz(temp));
            members.push_back(chef);
            temp &= (temp - 1U);
        }

        const unsigned m = static_cast<unsigned>(members.size());
        if (m == 1U) {
            const unsigned chef = members[0];
            const std::size_t s = state_index(mask, chef, n);
            wins[s * win_stride + chef] = 1.0L;
            dishes[s] = 0.0L;
            return;
        }

        std::vector<long double> a_win(static_cast<std::size_t>(m) * win_stride, 0.0L);
        std::vector<long double> b(static_cast<std::size_t>(m), 0.0L);
        std::vector<long double> a_dish(static_cast<std::size_t>(m), 0.0L);

        constexpr long double kEps = 1e-18L;

        for (unsigned t = 0U; t < m; ++t) {
            const unsigned current = members[t];
            const long double p = skills[current];
            b[t] = 1.0L - p;

            long double best_self_win = -1.0L;
            unsigned best_target_pos = 0U;
            unsigned best_distance = std::numeric_limits<unsigned>::max();

            for (unsigned u = 0U; u < m; ++u) {
                if (u == t) {
                    continue;
                }

                const unsigned target = members[u];
                const u32 next_mask = mask & ~(1U << target);
                const unsigned next_turn = next_alive_after(next_mask, current);
                const std::size_t succ_state = state_index(next_mask, next_turn, n);
                const long double self_win = wins[succ_state * win_stride + current];

                const unsigned distance = (u + m - t) % m;
                if (self_win > best_self_win + kEps ||
                    (std::fabsl(self_win - best_self_win) <= kEps &&
                     distance < best_distance)) {
                    best_self_win = self_win;
                    best_target_pos = u;
                    best_distance = distance;
                }
            }

            const unsigned target = members[best_target_pos];
            const u32 next_mask = mask & ~(1U << target);
            const unsigned next_turn = next_alive_after(next_mask, current);
            const std::size_t succ_state = state_index(next_mask, next_turn, n);

            long double* const a_row = &a_win[static_cast<std::size_t>(t) * win_stride];
            const long double* const succ_row = &wins[succ_state * win_stride];
            for (unsigned k = 0U; k < n; ++k) {
                a_row[k] = p * succ_row[k];
            }
            a_dish[t] = 1.0L + p * dishes[succ_state];
        }

        std::vector<long double> x0(win_stride, 0.0L);
        std::vector<long double> acc(win_stride, 0.0L);
        long double B = b[0];
        {
            const long double* row0 = &a_win[0];
            for (unsigned k = 0U; k < n; ++k) {
                acc[k] = row0[k];
            }
        }

        for (unsigned t = 1U; t < m; ++t) {
            const long double* row = &a_win[static_cast<std::size_t>(t) * win_stride];
            for (unsigned k = 0U; k < n; ++k) {
                acc[k] += B * row[k];
            }
            B *= b[t];
        }

        const long double den = 1.0L - B;
        if (std::fabsl(den) < 1e-24L) {
            throw std::runtime_error("Numerical instability: denominator close to zero.");
        }
        for (unsigned k = 0U; k < n; ++k) {
            x0[k] = acc[k] / den;
        }

        std::vector<long double> state_wins(static_cast<std::size_t>(m) * win_stride, 0.0L);
        long double* const row_start = &state_wins[0];
        for (unsigned k = 0U; k < n; ++k) {
            row_start[k] = x0[k];
        }

        long double* const row_last = &state_wins[static_cast<std::size_t>(m - 1U) * win_stride];
        const long double* const a_last = &a_win[static_cast<std::size_t>(m - 1U) * win_stride];
        for (unsigned k = 0U; k < n; ++k) {
            row_last[k] = a_last[k] + b[m - 1U] * x0[k];
        }
        for (int t = static_cast<int>(m) - 2; t >= 1; --t) {
            long double* const row = &state_wins[static_cast<std::size_t>(t) * win_stride];
            const long double* const a_row = &a_win[static_cast<std::size_t>(t) * win_stride];
            const long double* const next_row =
                &state_wins[static_cast<std::size_t>(t + 1) * win_stride];
            for (unsigned k = 0U; k < n; ++k) {
                row[k] = a_row[k] + b[static_cast<std::size_t>(t)] * next_row[k];
            }
        }

        long double B_d = b[0];
        long double A_d = a_dish[0];
        for (unsigned t = 1U; t < m; ++t) {
            A_d += B_d * a_dish[t];
            B_d *= b[t];
        }
        const long double den_d = 1.0L - B_d;
        if (std::fabsl(den_d) < 1e-24L) {
            throw std::runtime_error("Numerical instability: denominator close to zero.");
        }
        const long double d0 = A_d / den_d;

        std::vector<long double> state_dishes(static_cast<std::size_t>(m), 0.0L);
        state_dishes[0] = d0;
        state_dishes[m - 1U] = a_dish[m - 1U] + b[m - 1U] * d0;
        for (int t = static_cast<int>(m) - 2; t >= 1; --t) {
            state_dishes[static_cast<std::size_t>(t)] =
                a_dish[static_cast<std::size_t>(t)] +
                b[static_cast<std::size_t>(t)] * state_dishes[static_cast<std::size_t>(t + 1)];
        }

        for (unsigned t = 0U; t < m; ++t) {
            const unsigned turn = members[t];
            const std::size_t s = state_index(mask, turn, n);
            const long double* const from_row =
                &state_wins[static_cast<std::size_t>(t) * win_stride];
            long double* const to_row = &wins[s * win_stride];
            for (unsigned k = 0U; k < n; ++k) {
                to_row[k] = from_row[k];
            }
            dishes[s] = state_dishes[t];
        }
    };

    for (unsigned size = 1U; size <= n; ++size) {
        const auto& masks = masks_by_size[size];
        const unsigned threads =
            choose_thread_count(allow_multithreading, requested_threads, masks.size());
        if (threads == 1U) {
            for (const u32 mask : masks) {
                process_mask(mask);
            }
            continue;
        }

        std::atomic<std::size_t> cursor{0U};
        std::vector<std::thread> pool;
        pool.reserve(threads);

        for (unsigned t = 0U; t < threads; ++t) {
            pool.emplace_back([&]() {
                while (true) {
                    const std::size_t idx = cursor.fetch_add(1U, std::memory_order_relaxed);
                    if (idx >= masks.size()) {
                        break;
                    }
                    process_mask(masks[idx]);
                }
            });
        }
        for (auto& th : pool) {
            th.join();
        }
    }

    CompetitionResult result;
    const u32 full_mask = total_masks - 1U;
    const std::size_t start_state = state_index(full_mask, 0U, n);
    result.win_from_start.assign(win_stride, 0.0L);
    {
        const long double* row = &wins[start_state * win_stride];
        for (unsigned k = 0U; k < n; ++k) {
            result.win_from_start[k] = row[k];
        }
    }
    result.expected_dishes_from_start = dishes[start_state];
    return result;
}

bool almost_equal(const long double a, const long double b, const long double eps) {
    return std::fabsl(a - b) <= eps;
}

bool run_checkpoints(const bool allow_multithreading, const unsigned requested_threads) {
    {
        const std::vector<long double> skills = {0.25L, 0.5L, 1.0L};
        const CompetitionResult result =
            solve_competition(skills, allow_multithreading, requested_threads);
        const long double expected_w31 = 0.29375L;
        if (!almost_equal(result.win_from_start[0], expected_w31, 1e-14L)) {
            std::cerr << "Checkpoint failed: W3(1). Expected " << std::setprecision(15)
                      << static_cast<double>(expected_w31) << ", got "
                      << static_cast<double>(result.win_from_start[0]) << '\n';
            return false;
        }
    }

    {
        const std::vector<long double> skills = fibonacci_skills(7U);
        const CompetitionResult result =
            solve_competition(skills, allow_multithreading, requested_threads);
        const long double expected_w[7] = {
            0.08965042L, 0.20775702L, 0.15291406L, 0.14554098L,
            0.15905291L, 0.10261412L, 0.14247050L};
        const long double expected_e7 = 42.28176050L;

        for (unsigned i = 0U; i < 7U; ++i) {
            const long double got = result.win_from_start[i];
            if (!almost_equal(got, expected_w[i], 5e-9L)) {
                std::cerr << "Checkpoint failed: W7(" << (i + 1U) << "). Expected "
                          << std::fixed << std::setprecision(8)
                          << static_cast<double>(expected_w[i]) << ", got "
                          << static_cast<double>(got) << '\n';
                return false;
            }
        }
        if (!almost_equal(result.expected_dishes_from_start, expected_e7, 5e-9L)) {
            std::cerr << "Checkpoint failed: E(7). Expected " << std::fixed
                      << std::setprecision(8) << static_cast<double>(expected_e7)
                      << ", got " << static_cast<double>(result.expected_dishes_from_start)
                      << '\n';
            return false;
        }
    }

    return true;
}

} // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints &&
        !run_checkpoints(options.allow_multithreading, options.requested_threads)) {
        return 1;
    }

    const std::vector<long double> skills = fibonacci_skills(options.n);
    const CompetitionResult result =
        solve_competition(skills, options.allow_multithreading, options.requested_threads);

    std::cout << std::fixed << std::setprecision(8)
              << "E(" << options.n << ") = "
              << static_cast<double>(result.expected_dishes_from_start) << '\n';
    return 0;
}
