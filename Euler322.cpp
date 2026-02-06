#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

std::vector<int> to_base_digits(u64 x, const int base) {
    std::vector<int> digits;
    if (x == 0ULL) {
        digits.push_back(0);
        return digits;
    }

    while (x > 0ULL) {
        digits.push_back(static_cast<int>(x % static_cast<u64>(base)));
        x /= static_cast<u64>(base);
    }
    return digits;
}

u64 pow_u64(u64 base, int exponent) {
    u64 result = 1ULL;
    while (exponent-- > 0) {
        result *= base;
    }
    return result;
}

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
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

        u64 parsed = 0ULL;
        if (parse_u64_after_prefix(arg, "--threads=", parsed)) {
            if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
                std::cerr << "--threads is too large.\n";
                return false;
            }
            options.requested_threads = static_cast<unsigned>(parsed);
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

u64 count_no_carry_prime(const u64 n, const u64 limit, const int p) {
    // Counts x in [0, limit) such that adding n and x in base p has no carries.
    if (limit == 0ULL) {
        return 0ULL;
    }

    const u64 max_value = limit - 1ULL;
    std::vector<int> n_digits = to_base_digits(n, p);
    std::vector<int> x_digits = to_base_digits(max_value, p);

    const int len = std::max<int>(n_digits.size(), x_digits.size());
    n_digits.resize(len, 0);
    x_digits.resize(len, 0);

    u64 tight = 1ULL;
    u64 loose = 0ULL;

    for (int pos = len - 1; pos >= 0; --pos) {
        const int allowed_max = p - 1 - n_digits[pos];
        const int bound_digit = x_digits[pos];

        u64 next_tight = 0ULL;
        u64 next_loose = 0ULL;

        if (allowed_max >= 0) {
            next_loose += loose * static_cast<u64>(allowed_max + 1);

            for (int d = 0; d <= std::min(allowed_max, bound_digit); ++d) {
                if (d < bound_digit) {
                    next_loose += tight;
                } else {
                    next_tight += tight;
                }
            }
        }

        tight = next_tight;
        loose = next_loose;
    }

    return tight + loose;
}

void generate_low_values_rec(const std::vector<int>& allowed_digits,
                             const int pos,
                             const int split,
                             const u64 current,
                             const u64 place_value,
                             std::vector<u64>& out) {
    if (pos == split) {
        out.push_back(current);
        return;
    }

    for (int d = 0; d <= allowed_digits[pos]; ++d) {
        generate_low_values_rec(
            allowed_digits, pos + 1, split, current + static_cast<u64>(d) * place_value, place_value * 5ULL, out);
    }
}

std::vector<u64> generate_low_values(const std::vector<int>& allowed_digits, const int split) {
    std::vector<u64> out;
    generate_low_values_rec(allowed_digits, 0, split, 0ULL, 1ULL, out);
    return out;
}

std::vector<u64> generate_high_values(const std::vector<int>& allowed_digits,
                                      const int split,
                                      const u64 high_bound) {
    std::vector<int> bound_digits = to_base_digits(high_bound, 5);
    const int len = static_cast<int>(bound_digits.size());

    std::vector<int> high_allowed(len, 4);
    for (int i = 0; i < len; ++i) {
        const int original_pos = i + split;
        if (original_pos < static_cast<int>(allowed_digits.size())) {
            high_allowed[i] = allowed_digits[original_pos];
        }
    }

    std::vector<u64> pow5(len + 1, 1ULL);
    for (int i = 1; i <= len; ++i) {
        pow5[i] = pow5[i - 1] * 5ULL;
    }

    struct Node {
        int pos = 0;
        bool tight = true;
        u64 value = 0ULL;
    };

    std::vector<Node> stack;
    stack.push_back(Node());

    std::vector<u64> out;
    while (!stack.empty()) {
        const Node node = stack.back();
        stack.pop_back();

        if (node.pos == len) {
            out.push_back(node.value);
            continue;
        }

        const int digit_index = len - 1 - node.pos;
        const int limit_digit = node.tight ? bound_digits[digit_index] : 4;
        const int allowed_digit = high_allowed[digit_index];
        const int max_digit = std::min(limit_digit, allowed_digit);

        for (int d = max_digit; d >= 0; --d) {
            Node next;
            next.pos = node.pos + 1;
            next.tight = node.tight && (d == limit_digit);
            next.value = node.value + static_cast<u64>(d) * pow5[digit_index];
            stack.push_back(next);
        }
    }

    return out;
}

u64 count_overlap_coprime_10(const u64 n,
                             const u64 limit,
                             const bool allow_multithreading,
                             const unsigned requested_threads) {
    // Counts x in [0, limit) such that C(n+x, n) is coprime to 10.
    if (limit == 0ULL) {
        return 0ULL;
    }

    const u64 max_value = limit - 1ULL;
    std::vector<int> n_digits = to_base_digits(n, 5);
    std::vector<int> x_digits = to_base_digits(max_value, 5);
    const int len = std::max<int>(n_digits.size(), x_digits.size());
    n_digits.resize(len, 0);
    x_digits.resize(len, 0);

    std::vector<int> allowed_digits(len, 4);
    for (int i = 0; i < len; ++i) {
        allowed_digits[i] = 4 - n_digits[i];
    }

    int split = 0;
    while (split < len && x_digits[split] == allowed_digits[split]) {
        ++split;
    }

    const u64 step = pow_u64(5ULL, split);
    const u64 high_bound = max_value / step;

    const std::vector<u64> low_values = generate_low_values(allowed_digits, split);
    const std::vector<u64> high_values = generate_high_values(allowed_digits, split, high_bound);

    std::vector<u64> scaled_high_values;
    scaled_high_values.reserve(high_values.size());
    for (const u64 h : high_values) {
        scaled_high_values.push_back(h * step);
    }

    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, low_values.size());
    std::vector<std::thread> workers;
    workers.reserve(threads);
    std::vector<u64> partial(threads, 0ULL);

    auto worker = [&](const unsigned tid) {
        const std::size_t begin = (low_values.size() * tid) / threads;
        const std::size_t end = (low_values.size() * (tid + 1U)) / threads;

        u64 local = 0ULL;
        for (std::size_t i = begin; i < end; ++i) {
            const u64 low = low_values[i];
            for (const u64 scaled : scaled_high_values) {
                const u64 x = low + scaled;
                if ((x & n) == 0ULL) {
                    ++local;
                }
            }
        }
        partial[tid] = local;
    };

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back(worker, tid);
    }
    for (std::thread& worker_thread : workers) {
        worker_thread.join();
    }

    u64 total = 0ULL;
    for (const u64 value : partial) {
        total += value;
    }
    return total;
}

u64 solve_T(const u64 m,
            const u64 n,
            const bool allow_multithreading,
            const unsigned requested_threads) {
    const u64 interval = m - n;

    const u64 not_divisible_by_2 = count_no_carry_prime(n, interval, 2);
    const u64 not_divisible_by_5 = count_no_carry_prime(n, interval, 5);
    const u64 coprime_to_10 =
        count_overlap_coprime_10(n, interval, allow_multithreading, requested_threads);

    return interval - not_divisible_by_2 - not_divisible_by_5 + coprime_to_10;
}

u64 v_p_factorial(u64 x, const int p) {
    u64 total = 0ULL;
    while (x > 0ULL) {
        x /= static_cast<u64>(p);
        total += x;
    }
    return total;
}

u64 brute_force_T(const u64 m, const u64 n) {
    u64 count = 0ULL;
    for (u64 i = n; i < m; ++i) {
        const u64 vp2 =
            v_p_factorial(i, 2) - v_p_factorial(n, 2) - v_p_factorial(i - n, 2);
        const u64 vp5 =
            v_p_factorial(i, 5) - v_p_factorial(n, 5) - v_p_factorial(i - n, 5);
        if (vp2 > 0ULL && vp5 > 0ULL) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints(const Options& options) {
    struct Case {
        u64 m;
        u64 n;
        u64 expected;
    };

    const std::vector<Case> known_cases = {
        {1000000000ULL, 10000000ULL - 10ULL, 989697000ULL},
    };

    for (const Case& c : known_cases) {
        const u64 got = solve_T(c.m, c.n, options.allow_multithreading, options.requested_threads);
        if (got != c.expected) {
            std::cerr << "Checkpoint failed for T(" << c.m << ", " << c.n << "): got " << got
                      << ", expected " << c.expected << '\n';
            return false;
        }
    }

    const std::vector<std::pair<u64, u64>> brute_cases = {
        {100ULL, 10ULL},
        {160ULL, 47ULL},
        {220ULL, 90ULL},
    };

    for (const auto& c : brute_cases) {
        const u64 fast = solve_T(c.first, c.second, false, 1U);
        const u64 brute = brute_force_T(c.first, c.second);
        if (fast != brute) {
            std::cerr << "Brute checkpoint failed for T(" << c.first << ", " << c.second
                      << "): fast=" << fast << ", brute=" << brute << '\n';
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

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const u64 m = 1000000000000000000ULL;
    const u64 n = 1000000000000ULL - 10ULL;
    const u64 answer = solve_T(m, n, options.allow_multithreading, options.requested_threads);
    std::cout << answer << '\n';
    return 0;
}
