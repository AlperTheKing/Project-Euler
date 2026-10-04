#include <algorithm>
#include <array>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using Real = boost::multiprecision::cpp_dec_float_50;
using Basis = std::array<Real, 4>;
constexpr int TARGET = 100'000;
const Real EPS("1e-30");

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

int payoff(const int i, const int j) {
    if (i == j) return 0;
    const int winner = (i - j) % 2 != 0 ? std::min(i, j) : std::max(i, j);
    return (winner == i ? 1 : -1) * (2 * winner - 1);
}

struct Strategy {
    int first;
    std::vector<Real> probability;

    int last() const { return first + static_cast<int>(probability.size()) - 1; }
};

Real expected_payoff(const Strategy& strategy, const int i) {
    Real result = 0;
    for (int j = strategy.first; j <= strategy.last(); ++j) {
        result += payoff(i, j) * strategy.probability[j - strategy.first];
    }
    return result;
}

Strategy interval_strategy(const int first, const int last) {
    const int count = last - first + 1;
    require(count >= 3 && count % 2 != 0, "odd support size");
    std::vector<Basis> basis(count);
    basis[0][0] = 1;
    basis[1][0] = basis[1][1] = 1;
    basis[2][0] = basis[2][1] = basis[2][2] = 1;
    std::array<Basis, 2> difference{};
    difference[0][1] = difference[1][2] = 1;

    // For d_i = p_(i+1)-p_i, (2i+5)d_(i+2) = (2i-1)d_i - 4.
    for (int offset = 0; offset < count - 3; ++offset) {
        const int i = first + offset;
        Basis& d = difference[offset % 2];
        for (Real& value : d) value = value * (2 * i - 1) / (2 * i + 5);
        d[3] -= Real(4) / (2 * i + 5);
        for (int q = 0; q < 4; ++q) basis[offset + 3][q] = basis[offset + 2][q] + d[q];
    }

    std::array<Basis, 3> matrix{};
    for (int offset = 0; offset < count; ++offset) {
        for (int q = 0; q < 4; ++q) {
            matrix[0][q] += basis[offset][q];
            matrix[1][q] += Real(payoff(first, first + offset)) / (2 * first - 1) * basis[offset][q];
            matrix[2][q] += Real(payoff(first + 1, first + offset)) / (2 * first + 1) * basis[offset][q];
        }
    }
    matrix[0][3] = 1 - matrix[0][3];
    matrix[1][3] = -matrix[1][3];
    matrix[2][3] = -matrix[2][3];

    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 3; ++row) {
            if (abs(matrix[row][col]) > abs(matrix[pivot][col])) pivot = row;
        }
        std::swap(matrix[col], matrix[pivot]);
        require(matrix[col][col] != 0, "nonsingular boundary equations");
        const Real divisor = matrix[col][col];
        for (int q = col; q < 4; ++q) matrix[col][q] /= divisor;
        for (int row = 0; row < 3; ++row) {
            if (row == col) continue;
            const Real multiplier = matrix[row][col];
            for (int q = col; q < 4; ++q) matrix[row][q] -= multiplier * matrix[col][q];
        }
    }

    Strategy result{first, std::vector<Real>(count)};
    for (int offset = 0; offset < count; ++offset) {
        result.probability[offset] = basis[offset][3];
        for (int q = 0; q < 3; ++q) result.probability[offset] += basis[offset][q] * matrix[q][3];
    }
    return result;
}

bool lower_options_unprofitable(const Strategy& strategy) {
    for (int i = std::max(1, strategy.first - 2); i < strategy.first; ++i) {
        if (expected_payoff(strategy, i) > EPS) return false;
    }
    return true;
}

bool probabilities_nonnegative(const Strategy& strategy) {
    return std::all_of(strategy.probability.begin(), strategy.probability.end(),
        [](const Real& value) { return value >= -EPS; });
}

void check_strategy(const Strategy& strategy, const int n) {
    std::array<Real, 2> mass{}, weight{}, prefix_mass{}, prefix_weight{};
    for (int i = strategy.first; i <= strategy.last(); ++i) {
        const Real& p = strategy.probability[i - strategy.first];
        require(p >= -EPS, "nonnegative equilibrium probability");
        mass[i % 2] += p;
        weight[i % 2] += (2 * i - 1) * p;
    }
    require(abs(mass[0] + mass[1] - 1) < EPS, "normalized equilibrium");
    for (int i = strategy.first; i <= strategy.last(); ++i) {
        const int parity = i % 2;
        const Real& p = strategy.probability[i - strategy.first];
        const Real value = (2 * i - 1) * (prefix_mass[parity] + mass[1 - parity] - prefix_mass[1 - parity] + p)
            + prefix_weight[parity] - prefix_weight[1 - parity] - weight[parity];
        require(abs(value) < EPS, "zero payoff on equilibrium support");
        prefix_mass[parity] += p;
        prefix_weight[parity] += (2 * i - 1) * p;
    }
    require(lower_options_unprofitable(strategy), "unprofitable lower options");
    for (int i = strategy.last() + 1; i <= n; ++i) {
        require(expected_payoff(strategy, i) <= EPS, "unprofitable upper options");
    }
}

Strategy initial_strategy(const int n) {
    if (n <= 2) return {1, {Real(1)}};
    for (int count = 3; count <= n; count += 2) {
        for (const int last : {n, n - 1}) {
            if (last < count) continue;
            Strategy candidate = interval_strategy(last - count + 1, last);
            if (!probabilities_nonnegative(candidate) || !lower_options_unprofitable(candidate)) continue;
            if (last < n && expected_payoff(candidate, n) > EPS) continue;
            check_strategy(candidate, n);
            return candidate;
        }
    }
    throw std::runtime_error("No equilibrium support found");
}

Real advance(Strategy& strategy, const int n) {
    if (expected_payoff(strategy, n) <= EPS) return 0;
    int first = std::min(n - 2, strategy.first + (n - strategy.first) % 2);
    for (int attempt = 0; attempt <= n; ++attempt) {
        require(first >= 1 && first <= n - 2, "valid support endpoints");
        Strategy candidate = interval_strategy(first, n);
        if (!probabilities_nonnegative(candidate)) {
            first += 2;
        } else if (!lower_options_unprofitable(candidate)) {
            first -= 2;
        } else {
            check_strategy(candidate, n);
            strategy = std::move(candidate);
            return strategy.probability.back();
        }
    }
    throw std::runtime_error("Equilibrium support search did not converge");
}

struct Task {
    int first = 0;
    int last = 0;
    Real sum = 0;
    std::exception_ptr error;
};

void* sum_worker(void* argument) {
    Task& task = *static_cast<Task*>(argument);
    try {
        Strategy strategy = initial_strategy(task.first - 1);
        for (int n = task.first; n <= task.last; ++n) task.sum += advance(strategy, n);
    } catch (...) {
        task.error = std::current_exception();
    }
    return nullptr;
}

Real solve(const int n, unsigned thread_count) {
    if (n < 3) return 0;
    thread_count = std::min(thread_count, static_cast<unsigned>(n - 2));
    require(thread_count > 0, "positive thread count");
    std::vector<Task> tasks(thread_count);
    std::vector<pthread_t> threads(thread_count);
    unsigned created = 0;
    for (unsigned t = 0; t < thread_count; ++t) {
        tasks[t].first = 3 + static_cast<int>(static_cast<long long>(n - 2) * t / thread_count);
        tasks[t].last = 2 + static_cast<int>(static_cast<long long>(n - 2) * (t + 1) / thread_count);
        if (thread_count == 1) {
            sum_worker(&tasks[t]);
        } else {
            if (pthread_create(&threads[t], nullptr, sum_worker, &tasks[t]) != 0) break;
            ++created;
        }
    }
    bool joined = true;
    for (unsigned t = 0; t < created; ++t) joined = pthread_join(threads[t], nullptr) == 0 && joined;
    require(thread_count == 1 || created == thread_count, "pthread_create");
    require(joined, "pthread_join");
    Real sum = 0;
    for (const Task& task : tasks) {
        if (task.error) std::rethrow_exception(task.error);
        sum += task.sum;
    }
    return sum;
}

std::vector<Real> dense_equilibrium(const Strategy& strategy) {
    const int count = static_cast<int>(strategy.probability.size());
    std::vector<std::vector<Real>> matrix(count, std::vector<Real>(count + 1));
    for (int row = 0; row < count - 1; ++row) {
        for (int col = 0; col < count; ++col) matrix[row][col] = payoff(strategy.first + row, strategy.first + col);
    }
    std::fill(matrix.back().begin(), matrix.back().end(), Real(1));
    for (int col = 0; col < count; ++col) {
        int pivot = col;
        for (int row = col + 1; row < count; ++row) {
            if (abs(matrix[row][col]) > abs(matrix[pivot][col])) pivot = row;
        }
        std::swap(matrix[col], matrix[pivot]);
        require(matrix[col][col] != 0, "dense equilibrium pivot");
        const Real divisor = matrix[col][col];
        for (int q = col; q <= count; ++q) matrix[col][q] /= divisor;
        for (int row = 0; row < count; ++row) {
            if (row == col) continue;
            const Real multiplier = matrix[row][col];
            for (int q = col; q <= count; ++q) matrix[row][q] -= multiplier * matrix[col][q];
        }
    }
    std::vector<Real> result(count);
    for (int row = 0; row < count; ++row) result[row] = matrix[row][count];
    return result;
}

void run_tests(const unsigned thread_count) {
    Strategy strategy{1, {Real(1)}};
    Real sum = 0;
    const Real rounded_tolerance("5e-11");
    for (int n = 2; n <= 200; ++n) {
        const Real p = advance(strategy, n);
        if (n >= 3) sum += p;
        const std::vector<Real> independent = dense_equilibrium(strategy);
        for (std::size_t i = 0; i < independent.size(); ++i) {
            require(abs(independent[i] - strategy.probability[i]) < EPS, "dense payoff-matrix comparison");
        }
        for (int i = 1; i <= n; ++i) require(expected_payoff(strategy, i) <= EPS, "all pure counter-strategies");
        if (n == 3) require(abs(p - Real(1) / 9) < EPS, "P(3) = 1/9");
        if (n == 4) require(abs(p - Real(1) / 5) < EPS, "P(4) = 1/5");
        if (n == 8) require(p == 0, "unused final option for n=8");
        if (n == 10) {
            require(abs(p - Real("0.0479638009")) < rounded_tolerance, "P(10)");
            require(abs(sum - Real("1.1546112276")) < rounded_tolerance, "S(10)");
        }
        if (n == 100) require(abs(sum - Real("4.8779925686")) < rounded_tolerance, "S(100)");
    }
    require(abs(solve(1000, 1) - solve(1000, thread_count)) < EPS, "thread consistency");
    std::cout << "All checks passed.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const unsigned thread_count = std::max(1U, std::thread::hardware_concurrency());
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests(thread_count);
            return EXIT_SUCCESS;
        }
        if (argc != 1) throw std::invalid_argument("Usage: Euler1012 [--self-test]");
        std::cout << std::fixed << std::setprecision(10) << solve(TARGET, thread_count) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
