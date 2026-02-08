#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct Options {
    u64 limit = 100000000000000000ULL;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
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
            throw std::overflow_error("u64 overflow");
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
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
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
        if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
            throw std::overflow_error("unsigned overflow");
        }
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

u64 abs_i128(const i128 value) {
    return static_cast<u64>(value < 0 ? -value : value);
}

u64 isqrt_u128(const u128 x) {
    if (x == 0U) {
        return 0U;
    }

    long double approx = std::sqrt(static_cast<long double>(x));
    u64 r = static_cast<u64>(approx);

    while (static_cast<u128>(r + 1U) * static_cast<u128>(r + 1U) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > x) {
        --r;
    }

    return r;
}

u64 fourth_root_bound(const u64 n) {
    const u128 target = static_cast<u128>(4U) * static_cast<u128>(n);
    long double approx = std::pow(static_cast<long double>(target), 0.25L);
    u64 m = static_cast<u64>(approx) + 4U;

    auto pow4 = [](const u64 x) -> u128 {
        const u128 y = static_cast<u128>(x);
        return y * y * y * y;
    };

    while (m > 0U && pow4(m) > target) {
        --m;
    }
    while (pow4(m + 1U) <= target) {
        ++m;
    }

    return m + 2U;
}

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

u64 count_triplets_fast(const u64 limit, unsigned thread_count) {
    if (limit == 0U) {
        return 0U;
    }

    const u64 m_max = fourth_root_bound(limit);
    if (m_max < 2U) {
        return 0U;
    }

    thread_count = std::max(1U, std::min(thread_count, static_cast<unsigned>(m_max - 1U)));

    std::vector<u128> partial(thread_count, 0U);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    for (unsigned tid = 0U; tid < thread_count; ++tid) {
        workers.emplace_back([=, &partial]() {
            u128 local = 0U;

            for (u64 m = 2U + static_cast<u64>(tid); m <= m_max; m += static_cast<u64>(thread_count)) {
                const u64 mm = m * m;
                const u64 n_start = m / 2U + 1U;

                for (u64 n = n_start; n < m; ++n) {
                    if (std::gcd(m, n) != 1U) {
                        continue;
                    }

                    const u64 nn = n * n;

                    // Parameterization of rational points on u^2+v^2=5.
                    const i128 x_raw = static_cast<i128>(mm) - static_cast<i128>(4U) *
                                       static_cast<i128>(m) * static_cast<i128>(n) -
                                       static_cast<i128>(nn);
                    const u64 x_abs = abs_i128(x_raw);
                    const u64 y_abs = 2U * (mm + m * n - nn);
                    const u64 w = mm + nn;

                    u64 g = std::gcd(x_abs, y_abs);
                    g = std::gcd(g, w);

                    // g divisible by 5 corresponds to the 5x-lift duplicate representation.
                    if (g % 5U == 0U) {
                        continue;
                    }

                    u64 u = x_abs / g;
                    u64 v = y_abs / g;
                    if (u < v) {
                        std::swap(u, v);
                    }

                    // n > m/2 guarantees the canonical region W < V < U < 2W.
                    const u128 a_primitive = (static_cast<u128>(u) * static_cast<u128>(v)) / 2U;
                    if (a_primitive > static_cast<u128>(limit)) {
                        continue;
                    }

                    local += static_cast<u128>(limit) / a_primitive;
                }
            }

            partial[tid] = local;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u128 total = 0U;
    for (const u128 value : partial) {
        total += value;
    }

    return static_cast<u64>(total);
}

u64 count_triplets_bruteforce(const u64 limit) {
    u64 total = 0U;

    for (u64 a = 1U; a <= limit; ++a) {
        const u128 a2 = static_cast<u128>(a) * static_cast<u128>(a);
        const u128 b2_max = (static_cast<u128>(8U) * a2 - 1U) / 5U;
        const u64 b_max = isqrt_u128(b2_max);

        for (u64 b = a + 1U; b <= b_max; ++b) {
            const u128 b2 = static_cast<u128>(b) * static_cast<u128>(b);
            const u128 den = static_cast<u128>(5U) * b2 - static_cast<u128>(4U) * a2;
            if (den == 0U) {
                continue;
            }

            const u128 num = static_cast<u128>(4U) * a2 * b2;
            if (num % den != 0U) {
                continue;
            }

            const u128 c2 = num / den;
            const u64 c = isqrt_u128(c2);
            if (static_cast<u128>(c) * static_cast<u128>(c) != c2) {
                continue;
            }
            if (!(b < c && c < 2U * a)) {
                continue;
            }

            ++total;
        }
    }

    return total;
}

void run_checkpoints() {
    struct Checkpoint {
        u64 n;
        u64 expected;
    };

    constexpr std::array<Checkpoint, 2U> known{{
        {100U, 0U},
        {10000U, 106U},
    }};

    for (const Checkpoint& checkpoint : known) {
        const u64 got = count_triplets_fast(checkpoint.n, 1U);
        if (got != checkpoint.expected) {
            throw std::runtime_error("Checkpoint failed at N=" + std::to_string(checkpoint.n));
        }
    }

    constexpr u64 brute_n = 5000U;
    const u64 fast_small = count_triplets_fast(brute_n, 1U);
    const u64 brute_small = count_triplets_bruteforce(brute_n);
    if (fast_small != brute_small) {
        throw std::runtime_error("Fast/bruteforce mismatch at N=" + std::to_string(brute_n));
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }

        if (options.run_checkpoints) {
            run_checkpoints();
        }

        const unsigned threads = pick_thread_count(options.requested_threads);
        const u64 answer = count_triplets_fast(options.limit, threads);
        std::cout << answer << '\n';
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
