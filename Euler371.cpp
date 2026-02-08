#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int max_number = 999;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--max-number=", options.max_number)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_number >= 3 && (options.max_number % 2 == 1);
}

long double expected_steps_compressed(const int pair_count) {
    const int total_numbers = 2 * pair_count + 2;
    std::vector<long double> e0(static_cast<std::size_t>(pair_count + 1), 0.0L);
    std::vector<long double> e1(static_cast<std::size_t>(pair_count + 1), 0.0L);

    for (int k = pair_count; k >= 0; --k) {
        const long double den = static_cast<long double>(total_numbers - (k + 1));
        const long double grow = static_cast<long double>(2 * (pair_count - k));
        if (k == pair_count) {
            e1[static_cast<std::size_t>(k)] = static_cast<long double>(total_numbers) / den;
            e0[static_cast<std::size_t>(k)] =
                (static_cast<long double>(total_numbers) + e1[static_cast<std::size_t>(k)]) / den;
        } else {
            e1[static_cast<std::size_t>(k)] =
                (static_cast<long double>(total_numbers) +
                 grow * e1[static_cast<std::size_t>(k + 1)]) /
                den;
            e0[static_cast<std::size_t>(k)] =
                (static_cast<long double>(total_numbers) +
                 grow * e0[static_cast<std::size_t>(k + 1)] +
                 e1[static_cast<std::size_t>(k)]) /
                den;
        }
    }
    return e0[0];
}

std::vector<long double> solve_linear_system(std::vector<std::vector<long double>> matrix) {
    const int n = static_cast<int>(matrix.size());
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row) {
            if (std::fabsl(matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]) >
                std::fabsl(matrix[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(col)])) {
                pivot = row;
            }
        }
        std::swap(matrix[static_cast<std::size_t>(col)], matrix[static_cast<std::size_t>(pivot)]);
        const long double diag = matrix[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
        for (int j = col; j <= n; ++j) {
            matrix[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)] /= diag;
        }
        for (int row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const long double factor = matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
            if (factor == 0.0L) {
                continue;
            }
            for (int j = col; j <= n; ++j) {
                matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(j)] -=
                    factor * matrix[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)];
            }
        }
    }

    std::vector<long double> solution(static_cast<std::size_t>(n), 0.0L);
    for (int i = 0; i < n; ++i) {
        solution[static_cast<std::size_t>(i)] = matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(n)];
    }
    return solution;
}

long double expected_steps_explicit_states(const int pair_count) {
    const int total_numbers = 2 * pair_count + 2;
    std::vector<int> pow3(static_cast<std::size_t>(pair_count + 1), 1);
    for (int i = 1; i <= pair_count; ++i) {
        pow3[static_cast<std::size_t>(i)] = 3 * pow3[static_cast<std::size_t>(i - 1)];
    }

    const int states = 2 * pow3[static_cast<std::size_t>(pair_count)];
    std::vector<std::vector<long double>> matrix(
        static_cast<std::size_t>(states),
        std::vector<long double>(static_cast<std::size_t>(states + 1), 0.0L));

    const long double prob = 1.0L / static_cast<long double>(total_numbers);

    for (int idx = 0; idx < states; ++idx) {
        const int self_seen = idx & 1;
        const int code = idx >> 1;

        matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>(idx)] = 1.0L;
        matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>(states)] = 1.0L;

        matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>(idx)] -= prob;  // Draw 0.

        if (self_seen == 0) {
            const int next_idx = (code << 1) | 1;
            matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>(next_idx)] -= prob;
        }
        // If self_seen == 1 and we draw self, the game ends (absorbing win), so no matrix term.

        for (int p = 0; p < pair_count; ++p) {
            const int place = pow3[static_cast<std::size_t>(p)];
            const int state = (code / place) % 3;

            if (state == 0) {
                const int left_code = code + place;      // 0 -> 1
                const int right_code = code + 2 * place; // 0 -> 2
                matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>((left_code << 1) | self_seen)] -=
                    prob;
                matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>((right_code << 1) | self_seen)] -=
                    prob;
            } else if (state == 1) {
                // Draw left: no change. Draw right: win.
                matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>(idx)] -= prob;
            } else {
                // Draw right: no change. Draw left: win.
                matrix[static_cast<std::size_t>(idx)][static_cast<std::size_t>(idx)] -= prob;
            }
        }
    }

    const std::vector<long double> solution = solve_linear_system(std::move(matrix));
    const int initial_state = 0;  // No pair-side seen, self not seen.
    return solution[static_cast<std::size_t>(initial_state)];
}

long double solve(const int max_number) {
    const int pair_count = (max_number - 1) / 2;
    return expected_steps_compressed(pair_count);
}

bool run_checkpoints() {
    const long double expected_m1 = 38.0L / 9.0L;
    if (std::fabsl(expected_steps_compressed(1) - expected_m1) > 1e-15L) {
        std::cerr << "Checkpoint failed for pair_count=1 closed form" << '\n';
        return false;
    }

    const long double m2_fast = expected_steps_compressed(2);
    const long double m2_exact = expected_steps_explicit_states(2);
    if (std::fabsl(m2_fast - m2_exact) > 1e-12L) {
        std::cerr << "Checkpoint failed for explicit-state cross-check at pair_count=2" << '\n';
        return false;
    }

    const long double m3_fast = expected_steps_compressed(3);
    const long double m3_exact = expected_steps_explicit_states(3);
    if (std::fabsl(m3_fast - m3_exact) > 1e-12L) {
        std::cerr << "Checkpoint failed for explicit-state cross-check at pair_count=3" << '\n';
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
    std::cout << std::fixed << std::setprecision(8) << solve(options.max_number) << '\n';
    return 0;
}
