#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 5000000000000000ULL;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg.rfind("--threads=", 0U) == 0U) {
            const std::string tail = arg.substr(10);
            if (tail.empty()) {
                return false;
            }
            u64 parsed = 0;
            for (char c : tail) {
                if (c < '0' || c > '9') {
                    return false;
                }
                parsed = parsed * 10ULL + static_cast<u64>(c - '0');
                if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
                    return false;
                }
            }
            options.requested_threads = static_cast<unsigned>(parsed);
            continue;
        }
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 5ULL;
}

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) > r && static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > x) {
        --r;
    }
    return r;
}

inline u64 q_value(const u64 n) {
    return 2ULL * n * n + 2ULL * n + 1ULL;
}

u64 max_index_for_limit(const u64 limit) {
    if (limit <= 5ULL) {
        return 0ULL;
    }
    u64 m = (isqrt_u64(2ULL * limit - 1ULL) - 1ULL) / 2ULL;
    while (m > 0ULL && q_value(m) >= limit) {
        --m;
    }
    while (q_value(m + 1ULL) < limit) {
        ++m;
    }
    return m;
}

u64 solve(const u64 limit) {
    const u64 m = max_index_for_limit(limit);
    if (m == 0ULL) {
        return 0ULL;
    }

    std::vector<u64> values(static_cast<std::size_t>(m + 1ULL), 0ULL);
    std::vector<std::uint8_t> composite(static_cast<std::size_t>(m + 1ULL), 0U);
    for (u64 i = 1ULL; i <= m; ++i) {
        values[static_cast<std::size_t>(i)] = q_value(i);
    }

    u64 count = 0ULL;
    for (u64 i = 1ULL; i <= m; ++i) {
        if (composite[static_cast<std::size_t>(i)] == 0U) {
            ++count;
        }
        const u64 prime = values[static_cast<std::size_t>(i)];
        if (prime == 1ULL) {
            continue;
        }

        for (u64 j = i + prime; j <= m; j += prime) {
            composite[static_cast<std::size_t>(j)] = 1U;
            u64 v = values[static_cast<std::size_t>(j)];
            while (v % prime == 0ULL) {
                v /= prime;
            }
            values[static_cast<std::size_t>(j)] = v;
        }

        const u64 mod = (i + 1ULL) % prime;
        u64 start = (mod == 0ULL) ? prime : (prime - mod);
        if (start == i) {
            start += prime;
        }
        for (u64 j = start; j <= m; j += prime) {
            composite[static_cast<std::size_t>(j)] = 1U;
            u64 v = values[static_cast<std::size_t>(j)];
            while (v % prime == 0ULL) {
                v /= prime;
            }
            values[static_cast<std::size_t>(j)] = v;
        }
    }

    return count;
}

u64 brute_small(const u64 limit) {
    auto is_prime = [](u64 x) {
        if (x < 2ULL) {
            return false;
        }
        if ((x & 1ULL) == 0ULL) {
            return x == 2ULL;
        }
        for (u64 p = 3ULL; p * p <= x; p += 2ULL) {
            if (x % p == 0ULL) {
                return false;
            }
        }
        return true;
    };

    u64 count = 0ULL;
    for (u64 n = 1ULL;; ++n) {
        const u64 q = q_value(n);
        if (q >= limit) {
            break;
        }
        if (is_prime(q)) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(1000ULL) != brute_small(1000ULL)) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
        return false;
    }
    if (solve(1000000ULL) != brute_small(1000000ULL)) {
        std::cerr << "Checkpoint failed for limit=1,000,000" << '\n';
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
    (void)options.requested_threads;
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    std::cout << solve(options.limit) << '\n';
    return 0;
}
