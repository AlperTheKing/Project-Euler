#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u128 = unsigned __int128;

constexpr int MOD = 1'000'000'007;
constexpr int TARGET_N = 50;
constexpr int MAX_TERMS = (TARGET_N + 1) / 2;
constexpr int DELTA_LIMIT = 25;
constexpr int DELTA_SIZE = 2 * DELTA_LIMIT + 1;

struct Options {
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
    int target_n = TARGET_N;
};

struct Transition {
    std::array<std::vector<std::pair<int, u32>>, 10> by_residue;
};

int positive_mod(const int value, const int mod) {
    const int result = value % mod;
    return result < 0 ? result + mod : result;
}

u32 mod_add(const u32 a, const u32 b) {
    const u32 sum = a + b;
    return sum >= static_cast<u32>(MOD) ? sum - static_cast<u32>(MOD) : sum;
}

u32 mod_mul(const u32 a, const u32 b) {
    return static_cast<u32>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(MOD));
}

std::vector<u32> convolve(const std::vector<u32>& lhs, const std::vector<u32>& rhs) {
    std::vector<u32> result(lhs.size() + rhs.size() - 1, 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i] == 0) {
            continue;
        }
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            if (rhs[j] == 0) {
                continue;
            }
            const u32 add = mod_mul(lhs[i], rhs[j]);
            result[i + j] = mod_add(result[i + j], add);
        }
    }
    return result;
}

std::vector<u32> append_digit_set(const std::vector<u32>& poly, const int low_digit, const int high_digit) {
    std::vector<u32> next(poly.size() + static_cast<std::size_t>(high_digit), 0);
    for (std::size_t i = 0; i < poly.size(); ++i) {
        if (poly[i] == 0) {
            continue;
        }
        for (int digit = low_digit; digit <= high_digit; ++digit) {
            next[i + static_cast<std::size_t>(digit)] =
                mod_add(next[i + static_cast<std::size_t>(digit)], poly[i]);
        }
    }
    return next;
}

struct Solver {
    explicit Solver(const Options& options)
        : thread_count(resolve_thread_count(options)),
          binom(build_binom()),
          digit_sums(build_digit_sum_polynomials()),
          transition_index(build_transition_index()),
          transitions(build_transitions()) {}

    u32 solve(const int n) const {
        assert(0 <= n && n <= TARGET_N);

        const int stride_len = (MAX_TERMS + 1) * (MAX_TERMS + 1) * DELTA_SIZE;
        const int stride_a = (MAX_TERMS + 1) * DELTA_SIZE;
        const int stride_b = DELTA_SIZE;
        const auto index = [&](const int used, const int left_active, const int right_active, const int delta) {
            return used * stride_len + left_active * stride_a + right_active * stride_b +
                   (delta + DELTA_LIMIT);
        };

        std::vector<u32> dp(static_cast<std::size_t>(n + 1) * static_cast<std::size_t>(stride_len), 0);

        for (int left_terms = 1; left_terms <= MAX_TERMS; ++left_terms) {
            for (int right_terms = 1; right_terms + left_terms <= MAX_TERMS; ++right_terms) {
                const int base_length = left_terms + right_terms - 1;
                if (base_length > n) {
                    continue;
                }
                dp[static_cast<std::size_t>(index(base_length, left_terms, right_terms, 0))] = 1;
            }
        }

        for (int used = 0; used <= n; ++used) {
            for (int left_active = 0; left_active <= MAX_TERMS; ++left_active) {
                for (int right_active = 0; right_active + left_active <= MAX_TERMS; ++right_active) {
                    if (left_active == 0 && right_active == 0) {
                        continue;
                    }

                    const int next_length = used + left_active + right_active;
                    if (next_length > n) {
                        continue;
                    }

                    for (int delta = -DELTA_LIMIT; delta <= DELTA_LIMIT; ++delta) {
                        const u32 current =
                            dp[static_cast<std::size_t>(index(used, left_active, right_active, delta))];
                        if (current == 0) {
                            continue;
                        }

                        for (int left_next = 0; left_next <= left_active; ++left_next) {
                            const u32 choose_left = binom[left_active][left_next];
                            for (int right_next = 0; right_next <= right_active; ++right_next) {
                                const u32 choose_right = binom[right_active][right_next];
                                const u32 structure_weight =
                                    mod_mul(current, mod_mul(choose_left, choose_right));

                                const int tindex =
                                    transition_index[left_active][left_next][right_active][right_next];
                                const auto& group =
                                    transitions[static_cast<std::size_t>(tindex)]
                                        .by_residue[static_cast<std::size_t>(positive_mod(-delta, 10))];

                                for (const auto& [diff, ways] : group) {
                                    const int total = diff + delta;
                                    const int next_delta = total / 10;
                                    if (next_delta < -DELTA_LIMIT || next_delta > DELTA_LIMIT) {
                                        continue;
                                    }
                                    const std::size_t target =
                                        static_cast<std::size_t>(index(next_length,
                                                                       left_next,
                                                                       right_next,
                                                                       next_delta));
                                    dp[target] = mod_add(dp[target], mod_mul(structure_weight, ways));
                                }
                            }
                        }
                    }
                }
            }
        }

        u32 answer = 0;
        for (int used = 0; used <= n; ++used) {
            answer = mod_add(answer, dp[static_cast<std::size_t>(index(used, 0, 0, 0))]);
        }
        return answer;
    }

private:
    struct Task {
        int left_active = 0;
        int left_next = 0;
        int right_active = 0;
        int right_next = 0;
        int index = -1;
    };

    unsigned thread_count;
    std::array<std::array<u32, MAX_TERMS + 1>, MAX_TERMS + 1> binom{};
    std::array<std::array<std::vector<u32>, MAX_TERMS + 1>, MAX_TERMS + 1> digit_sums{};
    std::array<std::array<std::array<std::array<int, MAX_TERMS + 1>, MAX_TERMS + 1>, MAX_TERMS + 1>,
               MAX_TERMS + 1>
        transition_index{};
    std::vector<Transition> transitions;

    static unsigned resolve_thread_count(const Options& options) {
        if (!options.allow_multithreading) {
            return 1;
        }
        if (options.requested_threads != 0) {
            return std::max(1U, options.requested_threads);
        }
        const unsigned detected = std::thread::hardware_concurrency();
        return detected == 0 ? 1U : detected;
    }

    static std::array<std::array<u32, MAX_TERMS + 1>, MAX_TERMS + 1> build_binom() {
        std::array<std::array<u32, MAX_TERMS + 1>, MAX_TERMS + 1> result{};
        result[0][0] = 1;
        for (int n = 1; n <= MAX_TERMS; ++n) {
            result[n][0] = 1;
            result[n][n] = 1;
            for (int k = 1; k < n; ++k) {
                result[n][k] = result[n - 1][k - 1] + result[n - 1][k];
            }
        }
        return result;
    }

    static std::array<std::array<std::vector<u32>, MAX_TERMS + 1>, MAX_TERMS + 1>
    build_digit_sum_polynomials() {
        std::array<std::array<std::vector<u32>, MAX_TERMS + 1>, MAX_TERMS + 1> result{};

        std::array<std::vector<u32>, MAX_TERMS + 1> free_digits{};
        std::array<std::vector<u32>, MAX_TERMS + 1> leading_digits{};
        free_digits[0] = {1};
        leading_digits[0] = {1};

        for (int terms = 1; terms <= MAX_TERMS; ++terms) {
            free_digits[terms] = append_digit_set(free_digits[terms - 1], 0, 9);
            leading_digits[terms] = append_digit_set(leading_digits[terms - 1], 1, 9);
        }

        for (int active = 0; active <= MAX_TERMS; ++active) {
            for (int next = 0; next <= active; ++next) {
                result[active][next] = convolve(leading_digits[active - next], free_digits[next]);
            }
        }

        return result;
    }

    std::array<std::array<std::array<std::array<int, MAX_TERMS + 1>, MAX_TERMS + 1>, MAX_TERMS + 1>,
               MAX_TERMS + 1>
    build_transition_index() const {
        std::array<std::array<std::array<std::array<int, MAX_TERMS + 1>, MAX_TERMS + 1>, MAX_TERMS + 1>,
                   MAX_TERMS + 1>
            result{};
        int next_index = 0;
        for (int left_active = 0; left_active <= MAX_TERMS; ++left_active) {
            for (int left_next = 0; left_next <= left_active; ++left_next) {
                for (int right_active = 0; right_active + left_active <= MAX_TERMS; ++right_active) {
                    for (int right_next = 0; right_next <= right_active; ++right_next) {
                        result[left_active][left_next][right_active][right_next] = next_index++;
                    }
                }
            }
        }
        return result;
    }

    Transition build_single_transition(const Task& task) const {
        Transition transition;
        const auto& left_poly = digit_sums[task.left_active][task.left_next];
        const auto& right_poly = digit_sums[task.right_active][task.right_next];

        const int min_diff = -static_cast<int>(right_poly.size()) + 1;
        const int max_diff = static_cast<int>(left_poly.size()) - 1;
        std::vector<u32> diff_counts(static_cast<std::size_t>(max_diff - min_diff + 1), 0);

        for (int left_sum = 0; left_sum < static_cast<int>(left_poly.size()); ++left_sum) {
            if (left_poly[static_cast<std::size_t>(left_sum)] == 0) {
                continue;
            }
            for (int right_sum = 0; right_sum < static_cast<int>(right_poly.size()); ++right_sum) {
                if (right_poly[static_cast<std::size_t>(right_sum)] == 0) {
                    continue;
                }
                const int diff = left_sum - right_sum;
                const u32 ways =
                    mod_mul(left_poly[static_cast<std::size_t>(left_sum)],
                            right_poly[static_cast<std::size_t>(right_sum)]);
                diff_counts[static_cast<std::size_t>(diff - min_diff)] =
                    mod_add(diff_counts[static_cast<std::size_t>(diff - min_diff)], ways);
            }
        }

        for (int diff = min_diff; diff <= max_diff; ++diff) {
            const u32 ways = diff_counts[static_cast<std::size_t>(diff - min_diff)];
            if (ways == 0) {
                continue;
            }
            transition.by_residue[static_cast<std::size_t>(positive_mod(diff, 10))].push_back({diff, ways});
        }

        return transition;
    }

    std::vector<Transition> build_transitions() const {
        std::vector<Task> tasks;
        for (int left_active = 0; left_active <= MAX_TERMS; ++left_active) {
            for (int left_next = 0; left_next <= left_active; ++left_next) {
                for (int right_active = 0; right_active + left_active <= MAX_TERMS; ++right_active) {
                    for (int right_next = 0; right_next <= right_active; ++right_next) {
                        tasks.push_back({left_active,
                                         left_next,
                                         right_active,
                                         right_next,
                                         transition_index[left_active][left_next][right_active][right_next]});
                    }
                }
            }
        }

        std::vector<Transition> result(tasks.size());
        const unsigned workers =
            std::min<unsigned>(thread_count, static_cast<unsigned>(std::max<std::size_t>(1, tasks.size())));

        if (workers <= 1) {
            for (const Task& task : tasks) {
                result[static_cast<std::size_t>(task.index)] = build_single_transition(task);
            }
            return result;
        }

        std::vector<std::thread> pool;
        pool.reserve(workers);

        for (unsigned worker = 0; worker < workers; ++worker) {
            pool.emplace_back([&, worker]() {
                for (std::size_t i = worker; i < tasks.size(); i += workers) {
                    const Task& task = tasks[i];
                    result[static_cast<std::size_t>(task.index)] = build_single_transition(task);
                }
            });
        }

        for (std::thread& thread : pool) {
            thread.join();
        }

        return result;
    }
};

void usage() {
    std::cerr
        << "Usage:\n"
        << "  ./Euler990 [n] [--skip-checkpoints] [--single-thread] [--threads=N]\n";
}

Options parse_options(const int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            options.requested_threads = 1;
            continue;
        }
        if (arg.rfind("--threads=", 0) == 0) {
            options.requested_threads = static_cast<unsigned>(std::stoul(arg.substr(10)));
            continue;
        }
        if (!arg.empty() && arg[0] == '-') {
            usage();
            std::exit(EXIT_FAILURE);
        }

        options.target_n = std::stoi(arg);
    }

    if (options.target_n < 0 || options.target_n > TARGET_N) {
        std::cerr << "n must satisfy 0 <= n <= " << TARGET_N << ".\n";
        std::exit(EXIT_FAILURE);
    }

    return options;
}

void run_checkpoints(const Solver& solver) {
    assert(solver.solve(0) == 0);
    assert(solver.solve(1) == 0);
    assert(solver.solve(2) == 0);
    assert(solver.solve(3) == 9);
    assert(solver.solve(5) == 171);
    assert(solver.solve(7) == 4878);
}

}  // namespace

int main(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    const Solver solver(options);

    if (options.run_checkpoints) {
        run_checkpoints(solver);
    }

    std::cout << solver.solve(options.target_n) << '\n';
    return 0;
}
