#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct Options {
    u64 limit = 10000000000ULL;
    unsigned int threads = std::max(1U, std::thread::hardware_concurrency());
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_uint_after_prefix(const std::string& arg,
                             const std::string& prefix,
                             unsigned int& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned int>::max())) {
        return false;
    }

    value = static_cast<unsigned int>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_uint_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.threads == 0U) {
        options.threads = 1U;
    }

    if (options.limit > static_cast<u64>(std::numeric_limits<i64>::max())) {
        std::cerr << "--limit must be <= " << std::numeric_limits<i64>::max()
                  << " to keep intermediate signed states safe.\n";
        return false;
    }

    return true;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }

    std::string digits;
    while (value > 0U) {
        const int digit = static_cast<int>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

u64 isqrt_u128(const u128 n) {
    u128 lo = 0U;
    u128 hi = static_cast<u128>(std::numeric_limits<u64>::max());

    while (lo < hi) {
        const u128 mid = (lo + hi + 1U) >> 1;
        if (mid <= n / mid) {
            lo = mid;
        } else {
            hi = mid - 1U;
        }
    }

    return static_cast<u64>(lo);
}

u128 seed_first_area(const i64 p) {
    const u128 pp = static_cast<u128>(p);
    return 8U * pp * pp * pp + pp;
}

i64 max_seed(const u64 limit) {
    i64 lo = 0;
    i64 hi = 1;

    while (seed_first_area(hi) <= static_cast<u128>(limit)) {
        if (hi > std::numeric_limits<i64>::max() / 2) {
            break;
        }
        hi *= 2;
    }

    while (lo < hi) {
        const i64 mid = lo + (hi - lo + 1) / 2;
        if (seed_first_area(mid) <= static_cast<u128>(limit)) {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }

    return lo;
}

struct State {
    i64 p;
    i64 q;
    u64 area;
};

bool next_state(const State& current,
                const u64 limit,
                State& same_p,
                State& swap_branch,
                State& sign_branch) {
    const i128 p = static_cast<i128>(current.p);
    const i128 p2 = p * p;

    const i128 a = 8 * p2 + 1;
    const i128 b = 4 * p;
    const i128 c = 4 * p * (4 * p2 + 1);

    const i128 q1 = a * static_cast<i128>(current.q) + b * static_cast<i128>(current.area);
    const i128 s1 = c * static_cast<i128>(current.q) + a * static_cast<i128>(current.area);

    if (s1 <= 0 || static_cast<u128>(s1) > static_cast<u128>(limit)) {
        return false;
    }
    if (q1 <= 0 || q1 > static_cast<i128>(std::numeric_limits<i64>::max())) {
        return false;
    }

    const i64 nq = static_cast<i64>(q1);
    const u64 ns = static_cast<u64>(s1);

    same_p = State{current.p, nq, ns};
    swap_branch = State{nq, current.p, ns};
    sign_branch = State{nq, -current.p, ns};
    return true;
}

u128 sum_for_seed(const i64 seed, const u64 limit) {
    u128 sum = 0U;
    std::vector<State> stack;
    stack.reserve(256);
    stack.push_back(State{seed, 0, static_cast<u64>(seed)});

    while (!stack.empty()) {
        const State current = stack.back();
        stack.pop_back();

        State same_p{};
        State swap_branch{};
        State sign_branch{};
        if (!next_state(current, limit, same_p, swap_branch, sign_branch)) {
            continue;
        }

        sum += static_cast<u128>(same_p.area);

        stack.push_back(same_p);
        stack.push_back(swap_branch);
        stack.push_back(sign_branch);
    }

    return sum;
}

u128 solve(const u64 limit, const unsigned int requested_threads) {
    const i64 pmax = max_seed(limit);
    if (pmax <= 0) {
        return 0U;
    }

    const unsigned int workers = std::max(
        1U,
        std::min<unsigned int>(requested_threads, static_cast<unsigned int>(pmax)));

    std::atomic<i64> next_seed(1);
    std::vector<u128> partials(workers, 0U);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    for (unsigned int t = 0; t < workers; ++t) {
        threads.emplace_back([&, t]() {
            u128 local = 0U;
            while (true) {
                const i64 seed = next_seed.fetch_add(1, std::memory_order_relaxed);
                if (seed > pmax) {
                    break;
                }
                local += sum_for_seed(seed, limit);
            }
            partials[t] = local;
        });
    }

    for (std::thread& thread : threads) {
        thread.join();
    }

    u128 total = 0U;
    for (const u128 partial : partials) {
        total += partial;
    }
    return total;
}

u128 brute_sum_unordered(const u64 limit) {
    if (limit == 0U) {
        return 0U;
    }

    u128 total = 0U;
    const u64 max_b = 2U * limit;

    for (u64 b = 1U; b <= max_b; ++b) {
        const u64 max_c = max_b / b + 1U;
        if (max_c < b) {
            continue;
        }

        for (u64 c = b; c <= max_c; ++c) {
            const u128 m2 =
                static_cast<u128>(b) * b * c * c + static_cast<u128>(b) * b + static_cast<u128>(c) * c;
            const u64 m = isqrt_u128(m2);
            if (static_cast<u128>(m) * m != m2) {
                continue;
            }
            if ((m & 1U) != 0U) {
                continue;
            }

            const u64 area = m / 2U;
            if (area <= limit) {
                total += area;
            }
        }
    }

    return total;
}

void run_checkpoints(const Options& options) {
    {
        constexpr u64 b = 2U;
        constexpr u64 c = 8U;
        const u128 m2 = static_cast<u128>(b) * b * c * c +
                        static_cast<u128>(b) * b +
                        static_cast<u128>(c) * c;
        const u64 m = isqrt_u128(m2);
        if (static_cast<u128>(m) * m != m2 || m != 18U) {
            throw std::runtime_error("Checkpoint failed: sample triangle should have m=18.");
        }
        if (m / 2U != 9U) {
            throw std::runtime_error("Checkpoint failed: sample triangle area should be 9.");
        }
    }

    {
        constexpr u64 small_limit = 1000U;
        const u128 brute = brute_sum_unordered(small_limit);
        const u128 fast = solve(small_limit, 1U);
        if (brute != fast) {
            throw std::runtime_error(
                "Checkpoint failed: recurrence and brute-force mismatch at N=1000.");
        }
    }

    {
        constexpr u64 statement_limit = 1000000U;
        constexpr u64 statement_value = 18018206U;
        const u128 sample = solve(statement_limit, options.threads);
        if (sample != static_cast<u128>(statement_value)) {
            throw std::runtime_error("Checkpoint failed: S(10^6) mismatch.");
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    try {
        if (options.run_checkpoints) {
            run_checkpoints(options);
        }

        const u128 answer = solve(options.limit, options.threads);
        std::cout << "S(" << options.limit << ") = " << to_string_u128(answer)
                  << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }

    return 0;
}
