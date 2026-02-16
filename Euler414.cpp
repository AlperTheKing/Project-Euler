#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1000000000000000000ULL;  // last 18 digits

struct Options {
    int k_from = 2;
    int k_to = 300;
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
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--k-from=", options.k_from) ||
            parse_int_after_prefix(arg, "--k-to=", options.k_to)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.k_from >= 2 && options.k_to >= options.k_from;
}

u64 pow_u64(u64 base, int exp) {
    u64 r = 1;
    for (int i = 0; i < exp; ++i) {
        r *= base;
    }
    return r;
}

std::pair<int, int> next_state(int p, int q, int b) {
    const u64 b2 = static_cast<u64>(b) * static_cast<u64>(b);
    const u64 b4 = b2 * b2;
    u64 x = static_cast<u64>(p) * (b4 - 1ULL) +
            static_cast<u64>(q) * (b2 - 1ULL) * static_cast<u64>(b);

    std::array<int, 5> d{};
    for (int i = 0; i < 5; ++i) {
        d[static_cast<std::size_t>(i)] = static_cast<int>(x % static_cast<u64>(b));
        x /= static_cast<u64>(b);
    }
    std::sort(d.begin(), d.end());
    return {d[4] - d[0], d[3] - d[1]};
}

i64 ways(int p, int q, int b) {
    i64 t = 0;
    if (q == 0) {
        if (p == 0) {
            t = 1;                        // aaaaa
        } else {
            t += 120 / 24 * 2;            // aaaae, aeeee
            t += 120 / 6 * (p - 1);       // accce
        }
    } else if (p == q) {
        t += 120 / 2 / 2 * (p - 1);       // aacee
        t += 120 / 6 / 2 * 2;             // aaaee, aaeee
    } else {
        t += 120 / 2 * (q - 1);           // abcee
        t += 120 / 2 / 2;                 // abbee
        t += 120 / 6;                     // abeee
        t *= 2;

        if (p - 2 >= q) {
            t += 120 / 2 * (p - 1 - q) * 2;          // abbde, abdde
            t += 120 * (p - 1 - q) * (q - 1);        // abcde
        }
    }
    return static_cast<i64>(b - p) * t;
}

int dfs_depth(int p, int q, int b, std::vector<std::vector<int>>& memo) {
    int& cell = memo[static_cast<std::size_t>(p)][static_cast<std::size_t>(q)];
    if (cell != -1) {
        return cell;
    }
    const auto nxt = next_state(p, q, b);
    if (nxt.first == p && nxt.second == q) {
        cell = 0;
    } else {
        cell = dfs_depth(nxt.first, nxt.second, b, memo) + 1;
    }
    return cell;
}

u64 S_of_base(int b) {
    const u64 total = pow_u64(static_cast<u64>(b), 5);

    std::vector<std::vector<int>> memo(
        static_cast<std::size_t>(b),
        std::vector<int>(static_cast<std::size_t>(b), -1));

    u64 ans = 0ULL;
    for (int p = 1; p <= b - 1; ++p) {
        for (int q = 0; q <= p; ++q) {
            const i64 w = ways(p, q, b);
            const int d = dfs_depth(p, q, b, memo);
            ans = (ans + static_cast<u64>((u128)(w % static_cast<i64>(MOD)) * d % MOD)) % MOD;
        }
    }

    // +1 first-step contribution for all non-special numbers:
    // total numbers minus {0}, minus all-equal positive numbers (b-1), minus Kaprekar constant.
    ans = (ans + total - static_cast<u64>(b) - 1ULL) % MOD;
    return ans;
}

bool run_checkpoints() {
    if (S_of_base(15) != 5274369ULL) {
        std::cerr << "Checkpoint failed: S(15)\n";
        return false;
    }
    if (S_of_base(111) != 400668930299ULL) {
        std::cerr << "Checkpoint failed: S(111)\n";
        return false;
    }
    return true;
}

u64 solve_sum(int k_from, int k_to) {
    u64 ans = 0ULL;
    for (int k = k_from; k <= k_to; ++k) {
        const int b = 6 * k + 3;
        ans = (ans + S_of_base(b)) % MOD;
    }
    return ans;
}

struct SumWorkerArgs {
    int tid = 0;
    int thread_count = 1;
    int k_from = 2;
    int k_to = 300;
    u64 partial = 0ULL;
};

void* solve_sum_worker(void* raw) {
    auto* args = static_cast<SumWorkerArgs*>(raw);
    u64 local = 0ULL;
    for (int k = args->k_from + args->tid; k <= args->k_to; k += args->thread_count) {
        const int b = 6 * k + 3;
        local = (local + S_of_base(b)) % MOD;
    }
    args->partial = local;
    return nullptr;
}

u64 solve_sum_parallel(int k_from, int k_to) {
    const int total_k = k_to - k_from + 1;
    if (total_k <= 1) {
        return solve_sum(k_from, k_to);
    }

    long hw = sysconf(_SC_NPROCESSORS_ONLN);
    if (hw <= 1) {
        return solve_sum(k_from, k_to);
    }

    const int thread_count = std::min<int>(total_k, static_cast<int>(hw));
    std::vector<pthread_t> threads(static_cast<std::size_t>(thread_count));
    std::vector<SumWorkerArgs> args(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        args[static_cast<std::size_t>(t)].tid = t;
        args[static_cast<std::size_t>(t)].thread_count = thread_count;
        args[static_cast<std::size_t>(t)].k_from = k_from;
        args[static_cast<std::size_t>(t)].k_to = k_to;
        args[static_cast<std::size_t>(t)].partial = 0ULL;
    }

    bool create_failed = false;
    int started = 0;
    for (int t = 0; t < thread_count; ++t) {
        if (pthread_create(&threads[static_cast<std::size_t>(t)], nullptr,
                           solve_sum_worker, &args[static_cast<std::size_t>(t)]) != 0) {
            create_failed = true;
            break;
        }
        ++started;
    }

    for (int t = 0; t < started; ++t) {
        pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
    }

    if (create_failed) {
        return solve_sum(k_from, k_to);
    }

    u64 ans = 0ULL;
    for (const auto& a : args) {
        ans = (ans + a.partial) % MOD;
    }
    return ans;
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

    const u64 answer = solve_sum_parallel(options.k_from, options.k_to);
    std::cout << answer << '\n';
    return 0;
}
