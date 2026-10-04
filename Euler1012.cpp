#include <algorithm>
#include <array>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace {

using Real = boost::multiprecision::cpp_dec_float_50;
using u64 = std::uint64_t;
using Integer = boost::multiprecision::cpp_int;
constexpr u64 TARGET = 100'000;
constexpr u64 MAX_LIMIT = 10'000'000'000'000'000ULL;
constexpr unsigned SERIES_START = 64;
constexpr unsigned SERIES_ORDER = 8;
const Real EPS("1e-40");

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

u64 threshold(const u64 k) {
    return (4 * k * k * k + 5 * k) / 6 + 2;
}

Real probability(const u64 n, const u64 k) {
    Real product = 1;
    for (u64 i = 1; i <= k; ++i) product *= Real(2 * n - 4 * i - 1) / (2 * n - 4 * i + 1);
    return ((2 * k * k + 1) * product - 2 * (k * k - 1)) / (3 * (2 * k + 1));
}

Real block_sum(const u64 k, const u64 n, const std::vector<Real>& central) {
    const u64 first = threshold(k);
    const u64 last = std::min(n, threshold(k + 1) - 2);
    if (last < first) return 0;

    Real harmonic = (boost::math::digamma(Real(last) - Real("0.5"))
        - boost::math::digamma(Real(first) - Real("1.5"))) / 2;
    Real weighted = 0;
    for (u64 i = 1; i <= k; ++i) {
        const Real residue = 2 * central[i - 1] * (2 * k - 2 * i + 1) * central[k - i];
        weighted += residue * harmonic;
        if (i < k) {
            harmonic += Real(1) / (2 * first - 4 * i - 3) + Real(1) / (2 * first - 4 * i - 1)
                - Real(1) / (2 * last - 4 * i - 1) - Real(1) / (2 * last - 4 * i + 1);
        }
    }
    return (Real(last - first + 1) - Real(2 * k * k + 1) * weighted / 3) / (2 * k + 1);
}

struct Result {
    Real value = 0;
    Real error_bound = 0;
};

std::array<Real, SERIES_ORDER + 1> central_moments(const u64 k) {
    constexpr int coefficients[SERIES_ORDER + 1][SERIES_ORDER] = {
        {1}, {0}, {1, 2}, {1, 3, 2}, {3, 12, 10, -4},
        {6, 30, 40, 0, -16}, {15, 90, 150, 0, -104, 32},
        {36, 252, 525, 105, -560, -84, 272},
        {91, 728, 1792, 560, -2618, -616, 2248, -544}
    };
    const Real inverse = Real(1) / (k - 1);
    std::array<Real, SERIES_ORDER + 1> moments{};
    moments[0] = 1;
    for (unsigned j = 2; j <= SERIES_ORDER; ++j) {
        for (int degree = static_cast<int>(j) - 1; degree >= 0; --degree) {
            moments[j] = moments[j] * inverse + coefficients[j][degree];
        }
    }
    return moments;
}

Result series_block_sum(const u64 k, const u64 n) {
    const u64 first = threshold(k);
    const u64 last = std::min(n, threshold(k + 1) - 2);
    if (last < first) return {};
    const Real h = k - 1;
    const Real a = 2 * Real(first) - k - 2;
    const Real b = 2 * Real(last + 1) - k - 2;
    const Real ia = 1 / a, ib = 1 / b;
    const Real ua = h * ia, ub = h * ib;
    const Real ia2 = ia * ia, ib2 = ib * ib;
    const Real ia4 = ia2 * ia2, ib4 = ib2 * ib2;
    const Real ia6 = ia4 * ia2, ib6 = ib4 * ib2;
    const auto moments = central_moments(k);
    Real pa = 1, pb = 1, weighted = 0, bound = 0;
    for (unsigned j = 0; j <= SERIES_ORDER; ++j) {
        const Real integral = j == 0 ? Real(log(b / a) / 2) : Real((pa - pb) / (2 * j));
        const Real second = Real(j + 1) / 6;
        const Real fourth = Real((j + 1) * (j + 2) * (j + 3)) / 90;
        const Real sixth = Real((j + 1) * (j + 2) * (j + 3) * (j + 4) * (j + 5)) / 945;
        const Real remainder = sixth * (pa * ia6 - pb * ib6);
        const Real harmonic = integral + (pa * ia - pb * ib) / 2
            + second * (pa * ia2 - pb * ib2) - fourth * (pa * ia4 - pb * ib4) + remainder;
        weighted += moments[j] * harmonic;
        bound += abs(moments[j]) * remainder;
        pa *= ua;
        pb *= ub;
    }
    const Real rho = 3 * h / a;
    bound += Real(last - first + 1) / a * pow(rho, SERIES_ORDER + 1) / (1 - rho);
    const Real coefficient = Real(2 * k * k + 1) * (2 * k) / (3 * (2 * k + 1));
    return {Real(last - first + 1) / (2 * k + 1) - coefficient * weighted, coefficient * bound};
}

struct Task {
    u64 n = 0;
    unsigned index = 0;
    unsigned stride = 1;
    const std::vector<Real>* central = nullptr;
    std::vector<Result>* sums = nullptr;
    std::exception_ptr error;
};

void* sum_worker(void* argument) {
    Task& task = *static_cast<Task*>(argument);
    try {
        for (unsigned k = task.index + 1; k < task.sums->size(); k += task.stride) {
            (*task.sums)[k] = k < SERIES_START
                ? Result{block_sum(k, task.n, *task.central), 0}
                : series_block_sum(k, task.n);
        }
    } catch (...) {
        task.error = std::current_exception();
    }
    return nullptr;
}

Result evaluate(const u64 n, unsigned thread_count) {
    require(n <= MAX_LIMIT, "supported limit through 10^16");
    if (n < 3) return {};
    unsigned blocks = 1;
    while (threshold(blocks + 1) <= n) ++blocks;
    std::vector<Real> central(std::min(blocks, SERIES_START));
    std::vector<Result> sums(blocks + 1);
    central[0] = 1;
    for (unsigned j = 1; j < central.size(); ++j) central[j] = central[j - 1] * (2 * j - 1) / (2 * j);
    thread_count = std::min(thread_count, blocks);
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
    Result result;
    for (unsigned k = 1; k <= blocks; ++k) {
        result.value += sums[k].value;
        result.error_bound += sums[k].error_bound;
    }
    require(result.error_bound < Real("1e-24"), "series truncation bound");
    return result;
}

Real solve(const u64 n, const unsigned thread_count) {
    return evaluate(n, thread_count).value;
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

Integer threshold_numerator(const u64 k, const u64 n) {
    Integer numerator = 1, denominator = 1;
    for (u64 i = 1; i <= k; ++i) {
        numerator *= 2 * Integer(n) - 4 * i - 1;
        denominator *= 2 * Integer(n) - 4 * i + 1;
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
        if (threshold(k + 1) <= static_cast<u64>(n)) ++k;
        const int last = static_cast<u64>(n + 1) == threshold(k + 1) ? n - 1 : n;
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
    for (u64 block = 2; block <= 12; ++block) {
        std::vector<Real> central(block);
        central[0] = 1;
        for (u64 j = 1; j < block; ++j) central[j] = central[j - 1] * (2 * j - 1) / (2 * j);
        std::array<Real, SERIES_ORDER + 1> direct{};
        for (u64 j = 0; j < block; ++j) {
            const Real mass = central[j] * (2 * block - 2 * j - 1) * central[block - j - 1] / block;
            const Real z = 4 * Real(j) / (block - 1) - 1;
            Real power = 1;
            for (unsigned order = 0; order <= SERIES_ORDER; ++order) {
                direct[order] += mass * power;
                power *= z;
            }
        }
        const auto moments = central_moments(block);
        for (unsigned order = 0; order <= SERIES_ORDER; ++order) {
            require(abs(direct[order] - moments[order]) < EPS, "central moments versus partial-fraction residues");
        }
    }
    u64 largest_block = 1;
    while (threshold(largest_block + 1) <= MAX_LIMIT) ++largest_block;
    for (const u64 block : std::array<u64, 4>{64, 65, 1000, largest_block}) {
        std::vector<Real> central(block);
        central[0] = 1;
        for (u64 j = 1; j < block; ++j) central[j] = central[j - 1] * (2 * j - 1) / (2 * j);
        const u64 first = threshold(block), last = threshold(block + 1) - 2;
        for (const u64 end : {first, first + (last - first) / 2, last}) {
            const Result approximation = series_block_sum(block, end);
            require(abs(approximation.value - block_sum(block, end, central))
                < approximation.error_bound + Real("1e-35"), "bounded series versus exact harmonic block");
        }
    }
    require(abs(solve(1'000'000, 1) - solve(1'000'000, thread_count)) < EPS,
        "series thread consistency");
    std::cout << "All checks passed.\n";
}

unsigned logical_processor_count() {
#ifdef _WIN32
    return static_cast<unsigned>(GetActiveProcessorCount(ALL_PROCESSOR_GROUPS));
#else
    return std::thread::hardware_concurrency();
#endif
}

void print_table(const unsigned thread_count) {
    std::cout << "Hardware threads: " << thread_count << "; three runs per input.\n";
    std::cout << "| N | S(N) | Median seconds | Series error bound |\n"
              << "|---:|---:|---:|---:|\n";
    solve(TARGET, thread_count);
    u64 n = 1;
    for (unsigned exponent = 1; exponent <= 16; ++exponent) {
        n *= 10;
        std::array<double, 3> seconds{};
        Result result;
        for (unsigned repeat = 0; repeat < seconds.size(); ++repeat) {
            const auto start = std::chrono::steady_clock::now();
            const Result current = evaluate(n, thread_count);
            seconds[repeat] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (repeat != 0) require(abs(current.value - result.value) < EPS, "table repeat consistency");
            result = current;
        }
        std::sort(seconds.begin(), seconds.end());
        std::cout << "| 10^" << exponent << " | " << std::fixed << std::setprecision(10)
                  << result.value << " | " << std::setprecision(6) << seconds[1] << " | "
                  << std::scientific << std::setprecision(2) << result.error_bound << " |" << std::endl;
    }
}

}

int main(int argc, char* argv[]) {
    try {
        const unsigned thread_count = std::max(1U, logical_processor_count());
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests(thread_count);
            return EXIT_SUCCESS;
        }
        if (argc == 2 && std::string(argv[1]) == "--table") {
            print_table(thread_count);
            return EXIT_SUCCESS;
        }
        u64 n = TARGET;
        if (argc == 3 && std::string(argv[1]) == "--limit") {
            std::size_t consumed = 0;
            n = std::stoull(argv[2], &consumed);
            require(consumed == std::string(argv[2]).size(), "integer limit");
        } else if (argc != 1) {
            throw std::invalid_argument("Usage: Euler1012 [--self-test | --table | --limit N]");
        }
        std::cout << std::fixed << std::setprecision(10) << solve(n, thread_count) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
