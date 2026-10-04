#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
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
using Integer = boost::multiprecision::cpp_int;
constexpr int TARGET = 100'000;
const Real EPS("1e-40");

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

long long threshold(const int k) {
    return (4LL * k * k * k + 5LL * k) / 6 + 2;
}

Real probability(const int n, const int k) {
    Real product = 1;
    for (int i = 1; i <= k; ++i) product *= Real(2 * n - 4 * i - 1) / (2 * n - 4 * i + 1);
    return ((2 * k * k + 1) * product - 2 * (k * k - 1)) / (3 * (2 * k + 1));
}

Real block_sum(const int k, const int n, const std::vector<Real>& central) {
    const int first = static_cast<int>(threshold(k));
    const int last = static_cast<int>(std::min<long long>(n, threshold(k + 1) - 2));
    if (last < first) return 0;

    // H_m = sum_(j=1..m) 1/(2j-1); this difference sums 1/(2n-3) over the block.
    Real harmonic = (boost::math::digamma(Real(last) - Real("0.5"))
        - boost::math::digamma(Real(first) - Real("1.5"))) / 2;
    Real weighted = 0;
    for (int i = 1; i <= k; ++i) {
        const Real residue = 2 * central[i - 1] * (2 * k - 2 * i + 1) * central[k - i];
        weighted += residue * harmonic;
        if (i < k) {
            harmonic += Real(1) / (2 * first - 4 * i - 3) + Real(1) / (2 * first - 4 * i - 1)
                - Real(1) / (2 * last - 4 * i - 1) - Real(1) / (2 * last - 4 * i + 1);
        }
    }
    return (Real(last - first + 1) - Real(2 * k * k + 1) * weighted / 3) / (2 * k + 1);
}

struct Task {
    int n = 0;
    unsigned index = 0;
    unsigned stride = 1;
    const std::vector<Real>* central = nullptr;
    std::vector<Real>* sums = nullptr;
    std::exception_ptr error;
};

void* sum_worker(void* argument) {
    Task& task = *static_cast<Task*>(argument);
    try {
        for (unsigned k = task.index + 1; k < task.sums->size(); k += task.stride) {
            (*task.sums)[k] = block_sum(static_cast<int>(k), task.n, *task.central);
        }
    } catch (...) {
        task.error = std::current_exception();
    }
    return nullptr;
}

Real solve(const int n, unsigned thread_count) {
    if (n < 3) return 0;
    int blocks = 1;
    while (threshold(blocks + 1) <= n) ++blocks;
    std::vector<Real> central(blocks), sums(blocks + 1);
    central[0] = 1;
    for (int j = 1; j < blocks; ++j) central[j] = central[j - 1] * (2 * j - 1) / (2 * j);
    thread_count = std::min(thread_count, static_cast<unsigned>(blocks));
    require(thread_count > 0, "positive thread count");
    std::vector<Task> tasks(thread_count);
    std::vector<pthread_t> threads(thread_count);
    unsigned created = 0;
    for (unsigned t = 0; t < thread_count; ++t) {
        tasks[t] = {n, t, thread_count, &central, &sums, {}};
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
    for (const Task& task : tasks) if (task.error) std::rethrow_exception(task.error);
    Real result = 0;
    for (int k = 1; k <= blocks; ++k) result += sums[k];
    return result;
}

int payoff(const int i, const int j) {
    if (i == j) return 0;
    const int winner = (i - j) % 2 != 0 ? std::min(i, j) : std::max(i, j);
    return (winner == i ? 1 : -1) * (2 * winner - 1);
}

std::vector<Real> dense_equilibrium(const int first, const int last) {
    const int count = last - first + 1;
    std::vector<std::vector<Real>> matrix(count, std::vector<Real>(count + 1));
    for (int row = 0; row < count - 1; ++row) {
        for (int col = 0; col < count; ++col) matrix[row][col] = payoff(first + row, first + col);
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

Integer threshold_numerator(const int k, const int n) {
    Integer numerator = 1, denominator = 1;
    for (int i = 1; i <= k; ++i) {
        numerator *= 2 * n - 4 * i - 1;
        denominator *= 2 * n - 4 * i + 1;
    }
    return (2 * k * k + 1) * numerator - 2 * (k * k - 1) * denominator;
}

void run_tests(const unsigned thread_count) {
    int max_block = 1;
    while (threshold(max_block) <= TARGET) ++max_block;
    for (int k = 1; k <= max_block; ++k) {
        const int first = static_cast<int>(threshold(k));
        require(threshold_numerator(k, first - 1) < 0, "exact lower threshold sign");
        require(threshold_numerator(k, first) > 0, "exact upper threshold sign");
    }
    Real sum = 0;
    int k = 1;
    for (int n = 3; n <= 200; ++n) {
        if (threshold(k + 1) <= n) ++k;
        const int last = n + 1 == threshold(k + 1) ? n - 1 : n;
        const int first = last - 2 * k;
        const auto p = dense_equilibrium(first, last);
        Real mass = 0;
        for (const Real& value : p) {
            require(value > 0, "positive active probabilities");
            mass += value;
        }
        require(abs(mass - 1) < EPS, "normalized dense equilibrium");
        for (int i = 1; i <= n; ++i) {
            Real value = 0;
            for (int j = first; j <= last; ++j) value += payoff(i, j) * p[j - first];
            require(value <= EPS, "all pure counter-strategies");
        }
        const Real final_probability = last == n ? p.back() : Real(0);
        require(abs(final_probability - (last == n ? probability(n, k) : Real(0))) < EPS,
            "closed probability versus dense payoff matrix");
        sum += final_probability;
        if (n <= 10 || n == 21 || n == 47 || n == 88 || n == 100 || n == 200) {
            require(abs(sum - solve(n, 1)) < EPS, "harmonic block sum versus dense equilibria");
        }
    }
    require(abs(probability(3, 1) - Real(1) / 9) < EPS, "P(3) = 1/9");
    require(abs(probability(4, 1) - Real(1) / 5) < EPS, "P(4) = 1/5");
    require(abs(probability(10, 2) - Real("0.0479638009")) < Real("5e-11"), "P(10)");
    require(abs(solve(10, 1) - Real("1.1546112276")) < Real("5e-11"), "S(10)");
    require(abs(solve(100, 1) - Real("4.8779925686")) < Real("5e-11"), "S(100)");
    require(abs(solve(TARGET, 1) - solve(TARGET, thread_count)) < EPS, "thread consistency");
    std::cout << "All checks passed.\n";
}

}

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
