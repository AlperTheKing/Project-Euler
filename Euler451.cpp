#include <pthread.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    int limit = 20'000'000;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 3;
}

u32 mod_inverse(u32 a, u32 mod) {
    i64 x0 = 1;
    i64 y0 = 0;
    i64 x1 = 0;
    i64 y1 = 1;
    i64 aa = static_cast<i64>(a);
    i64 mm = static_cast<i64>(mod);

    while (mm != 0) {
        const i64 q = aa / mm;
        const i64 aa2 = aa - q * mm;
        aa = mm;
        mm = aa2;

        const i64 nx = x0 - q * x1;
        x0 = x1;
        x1 = nx;

        const i64 ny = y0 - q * y1;
        y0 = y1;
        y1 = ny;
    }

    if (aa != 1) {
        return 0U;
    }
    i64 v = x0 % static_cast<i64>(mod);
    if (v < 0) {
        v += static_cast<i64>(mod);
    }
    return static_cast<u32>(v);
}

std::vector<int> build_spf(const int limit) {
    std::vector<int> spf(static_cast<std::size_t>(limit + 1), 0);
    for (int i = 0; i <= limit; ++i) {
        spf[static_cast<std::size_t>(i)] = i;
    }
    for (int i = 2; i <= limit / i; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            if (spf[static_cast<std::size_t>(j)] == j) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return spf;
}

u32 compute_I(const int n, const std::vector<int>& spf) {
    std::array<u32, 256> sums{};
    const u32 target = static_cast<u32>(n - 1);
    sums[0] = target;
    int count = 1;
    u32 best = 1U;

    auto push = [&](u32 s) {
        sums[static_cast<std::size_t>(count++)] = s;
        if (s > best && s < target) {
            best = s;
        }
    };

    int nn = n;
    while (nn != 1) {
        const int p = spf[static_cast<std::size_t>(nn)];
        u32 q = 1U;
        do {
            nn /= p;
            q *= static_cast<u32>(p);
        } while (nn % p == 0);

        if (q == 2U) {
            continue;
        }

        const u32 m = static_cast<u32>(n / static_cast<int>(q));
        const u32 inv = mod_inverse(m % q, q);
        const u32 b = static_cast<u32>((static_cast<u64>(m) * inv) % static_cast<u64>(n));
        const u32 delta2 = static_cast<u32>((2ULL * b) % static_cast<u64>(n));
        u32 delta_half = 0U;
        u32 delta_half_plus_two = 0U;
        if ((q & 1U) == 0U && q >= 8U) {
            const u32 half = q / 2U;
            delta_half = static_cast<u32>((static_cast<u64>(half) * b) % static_cast<u64>(n));
            delta_half_plus_two =
                static_cast<u32>((static_cast<u64>(half + 2U) * b) % static_cast<u64>(n));
        }

        const int current_count = count;
        for (int j = 0; j < current_count; ++j) {
            const u32 base = sums[static_cast<std::size_t>(j)];

            u32 s = base + delta2;
            if (s >= static_cast<u32>(n)) {
                s -= static_cast<u32>(n);
            }
            push(s);

            if ((q & 1U) == 0U && q >= 8U) {
                u32 s1 = base + delta_half;
                if (s1 >= static_cast<u32>(n)) {
                    s1 -= static_cast<u32>(n);
                }
                push(s1);

                u32 s2 = base + delta_half_plus_two;
                if (s2 >= static_cast<u32>(n)) {
                    s2 -= static_cast<u32>(n);
                }
                push(s2);
            }
        }
    }
    return best;
}

u64 solve_range(const std::vector<int>& spf, const int begin, const int end) {
    u64 sum = 0ULL;
    for (int n = begin; n <= end; ++n) {
        sum += compute_I(n, spf);
    }
    return sum;
}

struct WorkerTask {
    const std::vector<int>* spf = nullptr;
    int begin = 0;
    int end = -1;
    u64 partial = 0ULL;
};

void* worker_entry(void* raw) {
    auto* task = static_cast<WorkerTask*>(raw);
    task->partial = solve_range(*task->spf, task->begin, task->end);
    return nullptr;
}

int detect_thread_count(const int work_items) {
    long cores = ::sysconf(_SC_NPROCESSORS_ONLN);
    int threads = (cores > 0) ? static_cast<int>(cores) : 4;
    if (threads < 1) {
        threads = 1;
    }
    if (threads > work_items) {
        threads = work_items;
    }
    return threads;
}

u64 solve(const int limit) {
    const std::vector<int> spf = build_spf(limit);
    const int begin = 3;
    if (limit < begin) {
        return 0ULL;
    }
    const int work_items = limit - begin + 1;
    const int threads = detect_thread_count(work_items);
    if (threads == 1) {
        return solve_range(spf, begin, limit);
    }

    std::vector<pthread_t> handles(static_cast<std::size_t>(threads));
    std::vector<WorkerTask> tasks(static_cast<std::size_t>(threads));
    std::vector<char> launched(static_cast<std::size_t>(threads), 0);

    const int chunk = (work_items + threads - 1) / threads;
    int start = begin;
    for (int t = 0; t < threads; ++t) {
        if (start > limit) {
            tasks[static_cast<std::size_t>(t)].partial = 0ULL;
            continue;
        }
        const int end = std::min(limit, start + chunk - 1);
        auto& task = tasks[static_cast<std::size_t>(t)];
        task.spf = &spf;
        task.begin = start;
        task.end = end;
        task.partial = 0ULL;

        const int rc = pthread_create(&handles[static_cast<std::size_t>(t)], nullptr, worker_entry, &task);
        if (rc == 0) {
            launched[static_cast<std::size_t>(t)] = 1;
        } else {
            task.partial = solve_range(spf, start, end);
        }
        start = end + 1;
    }

    u64 sum = 0ULL;
    for (int t = 0; t < threads; ++t) {
        if (launched[static_cast<std::size_t>(t)] != 0) {
            pthread_join(handles[static_cast<std::size_t>(t)], nullptr);
        }
        sum += tasks[static_cast<std::size_t>(t)].partial;
    }
    return sum;
}

u32 brute_I(const int n) {
    u32 best = 1U;
    for (int m = 1; m < n - 1; ++m) {
        if (std::gcd(m, n) != 1) {
            continue;
        }
        if ((static_cast<u64>(m) * static_cast<u64>(m)) % static_cast<u64>(n) == 1ULL) {
            best = static_cast<u32>(m);
        }
    }
    return best;
}

bool run_checkpoints() {
    const std::vector<int> spf = build_spf(2'000);
    if (compute_I(7, spf) != 1U) {
        std::cerr << "Checkpoint failed: I(7)=1" << '\n';
        return false;
    }
    if (compute_I(100, spf) != 51U) {
        std::cerr << "Checkpoint failed: I(100)=51" << '\n';
        return false;
    }
    for (int n = 3; n <= 600; ++n) {
        if (compute_I(n, spf) != brute_I(n)) {
            std::cerr << "Checkpoint failed: brute-force cross-check at n=" << n << '\n';
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
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    std::cout << solve(options.limit) << '\n';
    return 0;
}
