#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kTargetN = 1'234'567ULL;
constexpr long double kSqrt5 = 2.23606797749978969640917366873127623544L;
constexpr long double kPhi = (1.0L + kSqrt5) / 2.0L;
constexpr long double kLogPhi = 0.48121182505960344749775891342436842314L;
constexpr long double kTailDrop = 120.0L;

struct Options {
    u64 n = kTargetN;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
};

struct Window {
    u64 left = 1;
    u64 right = 1;
    u64 peak = 1;
    long double peak_log_term = 0.0L;
};

long double log_term(u64 n, u64 d) {
    return static_cast<long double>(n) * std::log(static_cast<long double>(d)) -
           static_cast<long double>(d) * kLogPhi;
}

unsigned resolve_thread_count(unsigned requested_threads) {
    if (requested_threads > 0) {
        return requested_threads;
    }
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 1;
    }
    return threads;
}

Window build_window(u64 n) {
    Window out;

    u64 approx = 1;
    if (n > 0) {
        const long double center = static_cast<long double>(n) / kLogPhi;
        if (center > static_cast<long double>(std::numeric_limits<u64>::max() - 10ULL)) {
            approx = std::numeric_limits<u64>::max() - 10ULL;
        } else {
            approx = static_cast<u64>(center);
            if (approx == 0) {
                approx = 1;
            }
        }
    }

    out.peak = approx;
    out.peak_log_term = log_term(n, out.peak);

    const u64 scan_lo = (approx > 6 ? approx - 6 : 1ULL);
    const u64 scan_hi = approx + 6;
    for (u64 d = scan_lo; d <= scan_hi; ++d) {
        const long double value = log_term(n, d);
        if (value > out.peak_log_term) {
            out.peak_log_term = value;
            out.peak = d;
        }
    }

    out.left = out.peak;
    while (out.left > 1) {
        const long double next_value = log_term(n, out.left - 1);
        if (out.peak_log_term - next_value > kTailDrop) {
            break;
        }
        --out.left;
    }

    out.right = out.peak;
    while (out.right < std::numeric_limits<u64>::max() - 1ULL) {
        const long double next_value = log_term(n, out.right + 1);
        if (out.peak_log_term - next_value > kTailDrop) {
            break;
        }
        ++out.right;
    }

    return out;
}

long double scaled_sum_range(u64 n,
                             u64 left,
                             u64 right,
                             long double peak_log_term,
                             unsigned threads) {
    if (left > right) {
        return 0.0L;
    }

    const u64 total_terms = right - left + 1ULL;
    if (threads <= 1 || total_terms < 50'000ULL) {
        long double total = 0.0L;
        for (u64 d = left; d <= right; ++d) {
            total += std::exp(log_term(n, d) - peak_log_term);
        }
        return total;
    }

    threads = std::min<unsigned>(threads, static_cast<unsigned>(total_terms));

    std::atomic<u64> next{left};
    std::vector<long double> partial(threads, 0.0L);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    constexpr u64 kChunk = 1024ULL;

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            long double local = 0.0L;
            while (true) {
                const u64 begin = next.fetch_add(kChunk, std::memory_order_relaxed);
                if (begin > right) {
                    break;
                }
                const u64 end = std::min(right, begin + kChunk - 1ULL);
                for (u64 d = begin; d <= end; ++d) {
                    local += std::exp(log_term(n, d) - peak_log_term);
                }
            }
            partial[tid] = local;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    long double total = 0.0L;
    for (long double value : partial) {
        total += value;
    }
    return total;
}

long double compute_log_series(u64 n, unsigned threads) {
    const Window window = build_window(n);

    const long double scaled =
        scaled_sum_range(n, window.left, window.right, window.peak_log_term, threads);

    return window.peak_log_term + std::log(scaled);
}

u64 compute_f(u64 n, unsigned threads) {
    const long double log_series = compute_log_series(n, threads);
    const long double raw = log_series / kLogPhi;

    // The closed form is F(n) = 3 + ceil(log_phi(sum_{d>=1} d^n / phi^d)).
    const long double rounded = std::ceil(raw - 1e-18L) + 3.0L;
    return static_cast<u64>(rounded);
}

bool run_validation_checkpoints(unsigned threads) {
    struct Checkpoint {
        u64 n;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {0ULL, 4ULL},
        {1ULL, 6ULL},
        {2ULL, 9ULL},
        {3ULL, 13ULL},
        {11ULL, 58ULL},
        {123ULL, 1173ULL},
    };

    for (const Checkpoint cp : checkpoints) {
        const u64 got = compute_f(cp.n, threads);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed: F(" << cp.n << ") expected " << cp.expected
                      << ", got " << got << "\n";
            return false;
        }
    }

    if (threads > 1) {
        const u64 thread_check_n = 200'000ULL;
        const u64 single = compute_f(thread_check_n, 1);
        const u64 multi = compute_f(thread_check_n, threads);
        if (single != multi) {
            std::cerr << "Thread consistency failed for n=" << thread_check_n << ": single="
                      << single << ", multi=" << multi << "\n";
            return false;
        }
    }

    std::cout << "Validation checkpoints passed.\n";
    return true;
}

bool parse_u64(const std::string& s, u64& out) {
    try {
        std::size_t pos = 0;
        const unsigned long long value = std::stoull(s, &pos, 10);
        if (pos != s.size()) {
            return false;
        }
        out = static_cast<u64>(value);
        return true;
    } catch (...) {
        return false;
    }
}

Options parse_options(int argc, char** argv) {
    Options options;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        if (arg.rfind("--n=", 0) == 0) {
            u64 value = 0;
            if (!parse_u64(arg.substr(4), value)) {
                std::cerr << "Invalid --n value: " << arg << "\n";
                std::exit(1);
            }
            options.n = value;
            continue;
        }

        if (arg.rfind("--threads=", 0) == 0) {
            u64 value = 0;
            if (!parse_u64(arg.substr(10), value) || value == 0ULL ||
                value > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
                std::cerr << "Invalid --threads value: " << arg << "\n";
                std::exit(1);
            }
            options.requested_threads = static_cast<unsigned>(value);
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        std::exit(1);
    }

    return options;
}

}  // namespace

int main(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    const unsigned threads = resolve_thread_count(options.requested_threads);

    if (options.run_checkpoints) {
        if (!run_validation_checkpoints(threads)) {
            return 1;
        }
    }

    const u64 answer = compute_f(options.n, threads);
    std::cout << answer << '\n';
    return 0;
}
