#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Options {
    int days = 10'000;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
    int initial_intervals = 4'096;
    int max_intervals = 1 << 20;
    long double relative_tolerance = 1e-12L;
};

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::int64_t parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const int digit = c - '0';
        if (parsed > (std::numeric_limits<std::int64_t>::max() - digit) / 10) {
            return false;
        }
        parsed = parsed * 10 + digit;
    }

    if (parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    int parsed = 0;
    if (!parse_int_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed < 0) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_long_double_after_prefix(const std::string& arg,
                                    const char* prefix,
                                    long double& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    char* end_ptr = nullptr;
    const long double parsed = std::strtold(tail.c_str(), &end_ptr);
    if (end_ptr == nullptr || *end_ptr != '\0' || !std::isfinite(parsed)) {
        return false;
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        int parsed_int = 0;
        if (parse_int_after_prefix(arg, "--days=", parsed_int) ||
            parse_int_after_prefix(arg, "--d=", parsed_int)) {
            options.days = parsed_int;
            continue;
        }
        if (parse_int_after_prefix(arg, "--initial-intervals=", parsed_int)) {
            options.initial_intervals = parsed_int;
            continue;
        }
        if (parse_int_after_prefix(arg, "--max-intervals=", parsed_int)) {
            options.max_intervals = parsed_int;
            continue;
        }

        unsigned parsed_unsigned = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        long double parsed_ld = 0.0L;
        if (parse_long_double_after_prefix(arg, "--rtol=", parsed_ld)) {
            options.relative_tolerance = parsed_ld;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.days < 2) {
        std::cerr << "--days must be >= 2.\n";
        return false;
    }
    if (options.initial_intervals < 4 || (options.initial_intervals % 2) != 0) {
        std::cerr << "--initial-intervals must be an even integer >= 4.\n";
        return false;
    }
    if (options.max_intervals < options.initial_intervals || (options.max_intervals % 2) != 0) {
        std::cerr << "--max-intervals must be an even integer >= --initial-intervals.\n";
        return false;
    }
    if (!(options.relative_tolerance > 0.0L) || !std::isfinite(options.relative_tolerance)) {
        std::cerr << "--rtol must be a positive finite number.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 20'000ULL) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }
    return std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

long double pow_int_long_double(long double base, int exponent) {
    if (exponent == 0) {
        return 1.0L;
    }

    int e = exponent;
    long double b = base;
    long double result = 1.0L;
    while (e > 0) {
        if ((e & 1) != 0) {
            result *= b;
        }
        b *= b;
        e >>= 1;
    }
    return result;
}

long double integrand_value(int days, long double x) {
    if (x <= 0.0L) {
        return 0.0L;
    }
    if (x >= 1.0L) {
        return 1.0L;
    }

    // Inclusion-exclusion + transfer matrix on the cycle gives:
    // E(D) / D = integral_0^1 (1 - Z_D(x)) / x dx,
    // Z_D(x) = lambda_1(x)^D + lambda_2(x)^D.
    const long double disc = std::sqrt((1.0L - x) * (1.0L + 3.0L * x));
    const long double lambda1 = (1.0L - x + disc) * 0.5L;
    const long double lambda2 = (1.0L - x - disc) * 0.5L;

    long double one_minus_lambda1_pow = 1.0L;
    if (lambda1 > 0.0L) {
        const long double log_term = static_cast<long double>(days) * std::log(lambda1);
        if (log_term > -1e-4L) {
            one_minus_lambda1_pow = -std::expm1(log_term);
        } else {
            one_minus_lambda1_pow = 1.0L - std::exp(log_term);
        }
    }

    const long double lambda2_pow = pow_int_long_double(lambda2, days);
    const long double numerator = one_minus_lambda1_pow - lambda2_pow;
    return numerator / x;
}

long double simpson_integral(int days,
                             int intervals,
                             bool allow_multithreading,
                             unsigned requested_threads) {
    const long double h = 1.0L / static_cast<long double>(intervals);
    const std::size_t workload = static_cast<std::size_t>(intervals - 1);
    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, workload);

    if (threads == 1) {
        long double sum = integrand_value(days, 0.0L) + integrand_value(days, 1.0L);
        for (int i = 1; i < intervals; ++i) {
            const long double x = static_cast<long double>(i) * h;
            const int weight = (i & 1) != 0 ? 4 : 2;
            sum += static_cast<long double>(weight) * integrand_value(days, x);
        }
        return static_cast<long double>(days) * h * sum / 3.0L;
    }

    std::vector<long double> partial(threads, 0.0L);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    const int block = (intervals - 1 + static_cast<int>(threads) - 1) / static_cast<int>(threads);
    for (unsigned t = 0; t < threads; ++t) {
        const int begin = 1 + static_cast<int>(t) * block;
        const int end = std::min(intervals, begin + block);
        workers.emplace_back([&, t, begin, end]() {
            long double local = 0.0L;
            for (int i = begin; i < end; ++i) {
                const long double x = static_cast<long double>(i) * h;
                const int weight = (i & 1) != 0 ? 4 : 2;
                local += static_cast<long double>(weight) * integrand_value(days, x);
            }
            partial[t] = local;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    long double sum = integrand_value(days, 0.0L) + integrand_value(days, 1.0L);
    for (const long double value : partial) {
        sum += value;
    }
    return static_cast<long double>(days) * h * sum / 3.0L;
}

long double compute_expectation(int days,
                                bool allow_multithreading,
                                unsigned requested_threads,
                                int initial_intervals,
                                int max_intervals,
                                long double relative_tolerance) {
    int intervals = initial_intervals;
    long double previous =
        simpson_integral(days, intervals, allow_multithreading, requested_threads);

    while (intervals < max_intervals) {
        intervals *= 2;
        const long double current =
            simpson_integral(days, intervals, allow_multithreading, requested_threads);
        const long double delta = std::fabs(current - previous);
        const long double scale = std::max(1.0L, std::fabs(current));
        if (delta <= relative_tolerance * scale) {
            return current;
        }
        previous = current;
    }

    return previous;
}

bool run_checkpoints(const Options& options) {
    const long double e2 = compute_expectation(
        2, options.allow_multithreading, options.requested_threads, 512, 32'768, 1e-13L);
    const long double e5 = compute_expectation(
        5, options.allow_multithreading, options.requested_threads, 512, 32'768, 1e-13L);
    const long double e365 = compute_expectation(
        365, options.allow_multithreading, options.requested_threads, 1'024, 65'536, 1e-12L);

    const long double expected_e2 = 1.0L;
    const long double expected_e5 = 31.0L / 6.0L;
    const long double expected_e365 = 1174.3501L;

    const bool ok2 = std::fabs(e2 - expected_e2) < 1e-11L;
    const bool ok5 = std::fabs(e5 - expected_e5) < 1e-11L;
    const bool ok365 = std::fabs(e365 - expected_e365) < 5e-5L;

    if (!ok2 || !ok5 || !ok365) {
        std::cerr << std::setprecision(18);
        std::cerr << "Checkpoint failed.\n";
        std::cerr << "E(2):   got " << e2 << ", expected " << expected_e2 << '\n';
        std::cerr << "E(5):   got " << e5 << ", expected " << expected_e5 << '\n';
        std::cerr << "E(365): got " << e365 << ", expected approx " << expected_e365 << '\n';
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

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const long double answer = compute_expectation(options.days,
                                                   options.allow_multithreading,
                                                   options.requested_threads,
                                                   options.initial_intervals,
                                                   options.max_intervals,
                                                   options.relative_tolerance);

    std::cout << std::fixed << std::setprecision(4) << answer << '\n';
    return 0;
}
