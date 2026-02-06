#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

namespace {

struct SolveResult {
    std::vector<std::vector<double>> f;      // P2 win prob at start of P2 turn
    std::vector<std::vector<int>> best_t;    // Optimal T at each (a,b)
    double start_probability = 0.0;          // P2 win prob from initial state (P1 starts)
};

unsigned compute_t_limit(int target) {
    unsigned t = 1;
    long long score = 1;
    while (score < target && t < 60) {
        ++t;
        score <<= 1;
    }
    return t + 4;
}

inline double p1_turn_value(const std::vector<std::vector<double>>& f, int a, int b) {
    if (a <= 0) return 0.0;
    if (b <= 0) return 1.0;
    const double prev_row = (a == 1) ? 0.0 : f[static_cast<std::size_t>(a - 1)][static_cast<std::size_t>(b)];
    return 0.5 * (f[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] + prev_row);
}

SolveResult solve_closed_form(int target) {
    SolveResult out;
    out.f.assign(static_cast<std::size_t>(target + 1), std::vector<double>(static_cast<std::size_t>(target + 1), 0.0));
    out.best_t.assign(static_cast<std::size_t>(target + 1), std::vector<int>(static_cast<std::size_t>(target + 1), 1));

    for (int a = 1; a <= target; ++a) {
        out.f[static_cast<std::size_t>(a)][0] = 1.0;
    }

    const unsigned t_limit = compute_t_limit(target);
    std::vector<double> q(static_cast<std::size_t>(t_limit + 1), 0.0);
    std::vector<int> score(static_cast<std::size_t>(t_limit + 1), 0);
    for (unsigned t = 1; t <= t_limit; ++t) {
        q[static_cast<std::size_t>(t)] = std::ldexp(1.0, -static_cast<int>(t));
        score[static_cast<std::size_t>(t)] = 1 << (t - 1);
    }

    for (int a = 1; a <= target; ++a) {
        for (int b = 1; b <= target; ++b) {
            const double prev_same_b = (a == 1) ? 0.0 : out.f[static_cast<std::size_t>(a - 1)][static_cast<std::size_t>(b)];

            double best_value = 0.0;
            int best_action = 1;

            for (unsigned t = 1; t <= t_limit; ++t) {
                const int gain = score[static_cast<std::size_t>(t)];
                const double success_value = (gain >= b)
                    ? 1.0
                    : p1_turn_value(out.f, a, b - gain);

                const double qt = q[static_cast<std::size_t>(t)];
                const double candidate = (2.0 * qt * success_value + (1.0 - qt) * prev_same_b) / (1.0 + qt);

                if (candidate > best_value) {
                    best_value = candidate;
                    best_action = static_cast<int>(t);
                }
            }

            out.f[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = best_value;
            out.best_t[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = best_action;
        }
    }

    out.start_probability = p1_turn_value(out.f, target, target);
    return out;
}

template <typename RowFn>
void parallel_rows(int row_begin, int row_end, unsigned threads, RowFn&& fn) {
    if (row_begin > row_end) return;
    if (threads <= 1) {
        for (int row = row_begin; row <= row_end; ++row) fn(row);
        return;
    }

    std::atomic<int> next_row{row_begin};
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&]() {
            while (true) {
                const int row = next_row.fetch_add(1, std::memory_order_relaxed);
                if (row > row_end) break;
                fn(row);
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }
}

// Independent Bellman value-iteration verifier.
double solve_value_iteration(int target, unsigned threads, double tolerance = 1e-14, int max_iterations = 5000) {
    if (target <= 0) return 0.0;

    const unsigned t_limit = compute_t_limit(target);
    std::vector<double> q(static_cast<std::size_t>(t_limit + 1), 0.0);
    std::vector<int> score(static_cast<std::size_t>(t_limit + 1), 0);
    for (unsigned t = 1; t <= t_limit; ++t) {
        q[static_cast<std::size_t>(t)] = std::ldexp(1.0, -static_cast<int>(t));
        score[static_cast<std::size_t>(t)] = 1 << (t - 1);
    }

    std::vector<std::vector<double>> p1(static_cast<std::size_t>(target + 1), std::vector<double>(static_cast<std::size_t>(target + 1), 0.0));
    std::vector<std::vector<double>> p2(static_cast<std::size_t>(target + 1), std::vector<double>(static_cast<std::size_t>(target + 1), 0.0));
    std::vector<std::vector<double>> next_p1(static_cast<std::size_t>(target + 1), std::vector<double>(static_cast<std::size_t>(target + 1), 0.0));
    std::vector<std::vector<double>> next_p2(static_cast<std::size_t>(target + 1), std::vector<double>(static_cast<std::size_t>(target + 1), 0.0));

    for (int a = 0; a <= target; ++a) {
        p1[static_cast<std::size_t>(a)][0] = 1.0;
        p2[static_cast<std::size_t>(a)][0] = 1.0;
        next_p1[static_cast<std::size_t>(a)][0] = 1.0;
        next_p2[static_cast<std::size_t>(a)][0] = 1.0;
    }

    std::vector<double> row_delta(static_cast<std::size_t>(target + 1), 0.0);

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        parallel_rows(1, target, threads, [&](int a) {
            for (int b = 1; b <= target; ++b) {
                const double p1_if_head = (a == 1) ? 0.0 : p2[static_cast<std::size_t>(a - 1)][static_cast<std::size_t>(b)];
                const double p1_if_tail = p2[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
                next_p1[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = 0.5 * (p1_if_head + p1_if_tail);
            }
        });

        parallel_rows(1, target, threads, [&](int a) {
            double local_max = 0.0;
            for (int b = 1; b <= target; ++b) {
                const double stay = next_p1[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
                double best_value = 0.0;

                for (unsigned t = 1; t <= t_limit; ++t) {
                    const int gain = score[static_cast<std::size_t>(t)];
                    const double success_value = (gain >= b)
                        ? 1.0
                        : next_p1[static_cast<std::size_t>(a)][static_cast<std::size_t>(b - gain)];
                    const double qt = q[static_cast<std::size_t>(t)];
                    const double candidate = qt * success_value + (1.0 - qt) * stay;
                    if (candidate > best_value) best_value = candidate;
                }

                next_p2[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = best_value;

                const double d1 = std::fabs(next_p1[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] -
                                            p1[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]);
                const double d2 = std::fabs(next_p2[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] -
                                            p2[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]);
                local_max = std::max(local_max, std::max(d1, d2));
            }
            row_delta[static_cast<std::size_t>(a)] = local_max;
        });

        double max_delta = 0.0;
        for (int a = 1; a <= target; ++a) {
            max_delta = std::max(max_delta, row_delta[static_cast<std::size_t>(a)]);
        }

        p1.swap(next_p1);
        p2.swap(next_p2);

        if (max_delta < tolerance) {
            break;
        }
    }

    return p1[static_cast<std::size_t>(target)][static_cast<std::size_t>(target)];
}

bool validate() {
    {
        const SolveResult tiny = solve_closed_form(1);
        const double expected = 1.0 / 3.0;
        if (std::fabs(tiny.start_probability - expected) > 1e-13) {
            std::cerr << "Validation failed for target=1: expected=" << expected
                      << ", got=" << tiny.start_probability << "\n";
            return false;
        }
    }

    const std::vector<int> deterministic_targets = {5, 10, 20};
    for (int target : deterministic_targets) {
        const double exact = solve_closed_form(target).start_probability;
        const double iter = solve_value_iteration(target, 1);
        if (std::fabs(exact - iter) > 2e-12) {
            std::cerr << "Validation failed for target=" << target
                      << ": exact=" << exact << ", value-iteration=" << iter << "\n";
            return false;
        }
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 2;
    const unsigned threads = std::min<unsigned>(hw, 8);
    if (threads > 1) {
        const int thread_check_target = 35;
        const double single = solve_value_iteration(thread_check_target, 1);
        const double multi = solve_value_iteration(thread_check_target, threads);
        if (std::fabs(single - multi) > 1e-12) {
            std::cerr << "Thread consistency failed for target=" << thread_check_target
                      << ": single=" << single << ", multi=" << multi << "\n";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!validate()) {
        return 1;
    }

    int target = 100;
    if (argc > 1) {
        target = std::max(1, std::atoi(argv[1]));
    }

    const SolveResult result = solve_closed_form(target);
    std::cout << std::fixed << std::setprecision(8) << result.start_probability << '\n';
    return 0;
}
