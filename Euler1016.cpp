#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
constexpr u64 MOD = 1'234'567'891;
constexpr u64 TARGET = 123'456'789;

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

constexpr u64 power(u64 base, u64 exponent) {
    u64 result = 1;
    while (exponent != 0) {
        if (exponent & 1) result = result * base % MOD;
        base = base * base % MOD;
        exponent >>= 1;
    }
    return result;
}

constexpr u64 INVERSE_15 = power(15, MOD - 2);
static_assert(15 * INVERSE_15 % MOD == 1);

struct Partial {
    u64 numerator = 1;
    u64 denominator = 1;
    u64 sum = 0;
};

Partial combine(const Partial& left, const Partial& right) {
    return {
        left.numerator * right.numerator % MOD,
        left.denominator * right.denominator % MOD,
        (left.sum * right.denominator + left.numerator * right.sum) % MOD
    };
}

struct Task {
    u64 first = 1;
    u64 last = 0;
    Partial result;
};

void* binomial_worker(void* argument) {
    auto& task = *static_cast<Task*>(argument);
    Partial value;
    // Common denominators make consecutive binomial sums composable without per-term division.
    for (u64 k = task.first; k <= task.last; ++k) {
        value.numerator = value.numerator * (4 * k - 2) % MOD;
        value.denominator = value.denominator * k % MOD;
        value.sum = (value.sum * k + value.numerator) % MOD;
    }
    task.result = value;
    return nullptr;
}

struct Binomials {
    u64 last;
    u64 sum;
};

Binomials central_binomials(const u64 limit, unsigned threads = 0) {
    if (limit == 0) return {1, 1};
    if (threads == 0) {
        threads = std::max(1U, std::thread::hardware_concurrency());
        threads = static_cast<unsigned>(std::min<u64>(threads, (limit + 999'999) / 1'000'000));
    }
    threads = static_cast<unsigned>(std::min<u64>(threads, limit));
    std::vector<Task> tasks(threads);
    std::vector<pthread_t> workers(threads);
    for (unsigned i = 0; i < threads; ++i) {
        tasks[i].first = limit * i / threads + 1;
        tasks[i].last = limit * (i + 1) / threads;
    }
    unsigned started = 1;
    for (; started < threads; ++started) {
        if (pthread_create(&workers[started], nullptr, binomial_worker, &tasks[started]) != 0) break;
    }
    binomial_worker(&tasks[0]);
    for (unsigned i = started; i < threads; ++i) binomial_worker(&tasks[i]);
    bool joined = true;
    for (unsigned i = 1; i < started; ++i) {
        if (pthread_join(workers[i], nullptr) != 0) joined = false;
    }
    require(joined, "pthread join");

    Partial total;
    for (const Task& task : tasks) total = combine(total, task.result);
    const u64 inverse = power(total.denominator, MOD - 2);
    require(total.denominator * inverse % MOD == 1, "invertible binomial denominator");
    return {total.numerator * inverse % MOD, (1 + total.sum * inverse % MOD) % MOD};
}

u64 solve(const u64 n, const unsigned threads = 0) {
    require(n < MOD, "limit below the modulus");
    if (n == 0) return 0;
    const auto binomial = central_binomials(n / 2, threads);
    const u64 square = n * n % MOD;
    const u64 cube = square * n % MOD;
    const u64 polynomial = (4 * cube + MOD - 12 * square % MOD + 50 * n
        + MOD - 63 + (n % 2 == 0 ? MOD - 1 : 1)) % MOD;
    const u64 coefficient = n % 2 != 0
        ? (16 * n % MOD) * ((2 * square + MOD - 8 * n % MOD + 11) % MOD) % MOD
        : (4 * (n + 1) % MOD) * ((4 * square + MOD - 18 * n % MOD + 35) % MOD) % MOD;
    const u64 numerator = (5 * power(2, n) % MOD * polynomial % MOD + 275
        + coefficient * binomial.last % MOD + MOD - 140 * binomial.sum % MOD) % MOD;
    return numerator * INVERSE_15 % MOD;
}

struct Enumeration {
    unsigned side;
    u64 count = 0;
    u64 traces = 0;

    void visit(const unsigned cell, const u64 used, const unsigned value,
               unsigned top, unsigned bottom, unsigned left, unsigned right,
               unsigned next_square, u64 trace) {
        const unsigned row = cell / side, column = cell % side;
        top = std::min(top, row);
        bottom = std::max(bottom, row);
        left = std::min(left, column);
        right = std::max(right, column);
        if (value == next_square * next_square) {
            if (bottom - top + 1 != next_square || right - left + 1 != next_square) return;
            ++next_square;
        }
        if (row == column) trace += value;
        if (value == side * side) {
            ++count;
            traces += trace;
            return;
        }
        const auto extend = [&](const unsigned next) {
            const u64 bit = u64{1} << next;
            if (!(used & bit)) visit(next, used | bit, value + 1, top, bottom, left, right,
                                    next_square, trace);
        };
        if (row != 0) extend(cell - side);
        if (row + 1 < side) extend(cell + side);
        if (column != 0) extend(cell - 1);
        if (column + 1 < side) extend(cell + 1);
    }
};

std::vector<u64> direct_trace_sums(const unsigned limit) {
    std::vector<u64> sums(limit + 1);
    if (limit == 0) return sums;
    sums[1] = 1;
    std::array<std::vector<u64>, 4> previous;
    for (auto& matrix : previous) matrix = {1};
    u64 count = 1;
    for (unsigned n = 2; n <= limit; ++n) {
        std::array<std::vector<u64>, 4> current;
        u64 trace = 0;
        for (unsigned corner = 0; corner < 4; ++corner) {
            current[corner].resize(n * n);
            const unsigned end_row = (corner / 2) * (n - 1);
            const unsigned end_column = (corner % 2) * (n - 1);
            for (unsigned dr = 0; dr < 2; ++dr) {
                const unsigned dc = corner == 0 || corner == 3 ? 1 - dr : dr;
                for (unsigned row = 0; row < n; ++row) {
                    for (unsigned column = 0; column < n; ++column) {
                        u64 value;
                        if (row >= dr && row < dr + n - 1
                            && column >= dc && column < dc + n - 1) {
                            value = previous[corner ^ 3][(row - dr) * (n - 1) + column - dc];
                        } else {
                            const unsigned distance = std::max(row, end_row) - std::min(row, end_row)
                                + std::max(column, end_column) - std::min(column, end_column);
                            value = count * (n * n - distance) % MOD;
                        }
                        u64& entry = current[corner][row * n + column];
                        entry = (entry + value) % MOD;
                    }
                }
            }
            for (unsigned i = 0; i < n; ++i) trace = (trace + current[corner][i * n + i]) % MOD;
        }
        sums[n] = (sums[n - 1] + trace) % MOD;
        previous = std::move(current);
        count = 2 * count % MOD;
    }
    return sums;
}

void run_tests() {
    for (u64 divisor = 2; divisor * divisor <= MOD; ++divisor) {
        require(MOD % divisor != 0, "prime modulus");
    }
    require(solve(0) == 0, "empty size range");
    require(solve(2) == 41, "given T(2)");
    require(solve(12) == 9'063'089, "given T(12)");
    require(solve(123) == 1'139'320'302, "given T(123)");

    u64 exhaustive_sum = 0;
    for (unsigned n = 1; n <= 6; ++n) {
        Enumeration enumeration{n};
        for (unsigned cell = 0; cell < n * n; ++cell) {
            enumeration.visit(cell, u64{1} << cell, 1, cell / n, cell / n, cell % n, cell % n, 1, 0);
        }
        require(enumeration.count == (n == 1 ? 1 : u64{1} << (n + 1)), "exhaustive matrix count");
        exhaustive_sum += enumeration.traces;
        require(solve(n, 1) == exhaustive_sum % MOD, "exhaustive Hamiltonian path traces");
    }

    const auto direct = direct_trace_sums(128);
    for (unsigned n = 0; n < direct.size(); ++n) {
        require(solve(n, 1) == direct[n], "direct matrix aggregation");
    }

    constexpr unsigned CHECK_LIMIT = 512;
    std::vector<u64> row(2 * CHECK_LIMIT + 1);
    row[0] = 1;
    u64 sum = 1;
    for (unsigned n = 1; n <= 2 * CHECK_LIMIT; ++n) {
        for (unsigned k = n; k > 0; --k) row[k] = (row[k] + row[k - 1]) % MOD;
        if (n % 2 != 0) continue;
        sum = (sum + row[n / 2]) % MOD;
        const auto actual = central_binomials(n / 2, 1);
        require(actual.last == row[n / 2] && actual.sum == sum, "Pascal triangle binomial sums");
    }
    const auto empty = central_binomials(0, 1);
    require(empty.last == 1 && empty.sum == 1, "empty binomial product");
    for (const u64 n : {1ULL, 2ULL, 7ULL, 31ULL, 128ULL, 511ULL, 1024ULL, 100'003ULL}) {
        const auto serial = central_binomials(n, 1);
        for (const unsigned threads : {2U, 3U, 7U}) {
            const auto parallel = central_binomials(n, threads);
            require(serial.last == parallel.last && serial.sum == parallel.sum,
                    "serial and pthread agreement");
        }
    }
    std::cout << "All checks passed; supplied examples, exhaustive matrices through size 6, "
                 "direct aggregation through size 128, Pascal binomial sums, and pthread agreement.\n";
}

u64 parse_integer(const std::string& text, const u64 maximum) {
    u64 value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || value > maximum) {
        throw std::invalid_argument("Expected an integer from 0 to " + std::to_string(maximum));
    }
    return value;
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") run_tests();
        else {
            if (argc > 3) throw std::invalid_argument("Usage: Euler1016 [n [threads] | --self-test]");
            const u64 n = argc >= 2 ? parse_integer(argv[1], MOD - 1) : TARGET;
            const unsigned threads = argc == 3 ? static_cast<unsigned>(parse_integer(argv[2], 256)) : 0;
            std::cout << solve(n, threads) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
