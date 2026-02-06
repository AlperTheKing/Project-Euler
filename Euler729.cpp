#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using i64 = std::int64_t;

struct Options {
    int p = 25;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct KahanSum {
    long double sum = 0.0L;
    long double c = 0.0L;

    void add(const long double x) {
        const long double y = x - c;
        const long double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
};

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::uint64_t parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
        return false;
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
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--p=", parsed_u32)) {
            if (parsed_u32 > static_cast<u32>(std::numeric_limits<int>::max())) {
                std::cerr << "--p is too large.\n";
                return false;
            }
            options.p = static_cast<int>(parsed_u32);
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.p < 2) {
        std::cerr << "--p must be at least 2.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload) {
    if (!allow_multithreading || workload < 2048ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

long double apply_inverse_word(long double y, const u32 bits, const int m) {
    for (int i = 0; i < m; ++i) {
        const long double d = std::sqrt(y * y + 4.0L);
        if (((bits >> i) & 1U) != 0U) {
            y = 0.5L * (y + d);
        } else {
            y = 0.5L * (y - d);
        }
    }
    return y;
}

long double solve_fixed_point(const u32 bits, const int m) {
    long double y = 0.0L;

    for (int it = 0; it < 80; ++it) {
        const long double ny = apply_inverse_word(y, bits, m);
        if (std::fabsl(ny - y) < 1e-28L) {
            return ny;
        }
        y = ny;
    }

    for (int it = 0; it < 16; ++it) {
        long double v = y;
        long double deriv = 1.0L;
        for (int i = 0; i < m; ++i) {
            const long double d = std::sqrt(v * v + 4.0L);
            const int s = (((bits >> i) & 1U) != 0U) ? 1 : -1;
            deriv *= 0.5L * (1.0L + static_cast<long double>(s) * (v / d));
            v = 0.5L * (v + static_cast<long double>(s) * d);
        }

        const long double phi = v - y;
        const long double dphi = deriv - 1.0L;
        if (std::fabsl(dphi) < 1e-30L) {
            break;
        }
        const long double step = phi / dphi;
        y -= step;
        if (std::fabsl(step) < 1e-28L) {
            break;
        }
    }

    return y;
}

long double orbit_range_from_inverse_word(const long double x0, const u32 bits, const int m) {
    long double mn = x0;
    long double mx = x0;
    long double y = x0;

    // The fixed point x0 of g_{s_{m-1}} ... g_{s_0} generates the same orbit
    // by repeatedly applying g_{s_0}, g_{s_1}, ..., g_{s_{m-2}}.
    for (int i = 0; i < m - 1; ++i) {
        const long double d = std::sqrt(y * y + 4.0L);
        if (((bits >> i) & 1U) != 0U) {
            y = 0.5L * (y + d);
        } else {
            y = 0.5L * (y - d);
        }
        mn = std::min(mn, y);
        mx = std::max(mx, y);
    }

    return mx - mn;
}

int mobius(const int n) {
    int x = n;
    int prime_count = 0;
    for (int p = 2; p * p <= x; ++p) {
        if (x % p != 0) {
            continue;
        }
        int exp = 0;
        while (x % p == 0) {
            x /= p;
            ++exp;
        }
        if (exp > 1) {
            return 0;
        }
        ++prime_count;
    }
    if (x > 1) {
        ++prime_count;
    }
    return (prime_count % 2 == 0) ? 1 : -1;
}

i64 expected_primitive_necklaces(const int n) {
    i64 sum = 0;
    for (int d = 1; d <= n; ++d) {
        if (n % d != 0) {
            continue;
        }
        const int mu = mobius(d);
        if (mu == 0) {
            continue;
        }
        sum += static_cast<i64>(mu) * (1LL << (n / d));
    }
    return sum / n;
}

void generate_primitive_necklaces(const int n, std::vector<u32>& out) {
    std::vector<int> a(static_cast<std::size_t>(n) + 1ULL, 0);
    out.clear();
    out.reserve(static_cast<std::size_t>(expected_primitive_necklaces(n)));

    std::function<void(int, int)> gen = [&](const int t, const int p) {
        if (t > n) {
            if (n % p == 0 && p == n) {
                u32 bits = 0U;
                for (int i = 1; i <= n; ++i) {
                    bits |= static_cast<u32>(a[i]) << static_cast<u32>(i - 1);
                }
                out.push_back(bits);
            }
            return;
        }

        a[t] = a[t - p];
        gen(t + 1, p);

        for (int j = a[t - p] + 1; j < 2; ++j) {
            a[t] = j;
            gen(t + 1, t);
        }
    };

    gen(1, 1);
}

long double solve_for_p(const int p,
                        const bool allow_multithreading,
                        const unsigned requested_threads,
                        const bool validate_counts) {
    KahanSum total;

    for (int m = 2; m <= p; ++m) {
        std::vector<u32> necklaces;
        generate_primitive_necklaces(m, necklaces);

        if (validate_counts) {
            const i64 expected = expected_primitive_necklaces(m);
            if (static_cast<i64>(necklaces.size()) != expected) {
                std::cerr << "Necklace count mismatch at m=" << m << ": expected " << expected
                          << ", got " << necklaces.size() << "\n";
                std::exit(1);
            }
        }

        const unsigned threads =
            choose_thread_count(allow_multithreading, requested_threads, necklaces.size());

        if (threads == 1U) {
            for (const u32 bits : necklaces) {
                const long double x0 = solve_fixed_point(bits, m);
                const long double r = orbit_range_from_inverse_word(x0, bits, m);
                total.add(static_cast<long double>(m) * r);
            }
            continue;
        }

        std::vector<KahanSum> partials(threads);
        auto worker = [&](const unsigned tid) {
            const std::size_t begin = (necklaces.size() * tid) / threads;
            const std::size_t end = (necklaces.size() * (tid + 1U)) / threads;
            for (std::size_t idx = begin; idx < end; ++idx) {
                const u32 bits = necklaces[idx];
                const long double x0 = solve_fixed_point(bits, m);
                const long double r = orbit_range_from_inverse_word(x0, bits, m);
                partials[tid].add(static_cast<long double>(m) * r);
            }
        };

        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (unsigned t = 0U; t < threads; ++t) {
            pool.emplace_back(worker, t);
        }
        for (auto& th : pool) {
            th.join();
        }

        for (const auto& ps : partials) {
            total.add(ps.sum);
            total.add(ps.c);
        }
    }

    return total.sum;
}

bool run_checkpoints(const bool allow_multithreading, const unsigned requested_threads) {
    struct Checkpoint {
        int p;
        long double expected;
    };

    const std::vector<Checkpoint> checks = {
        {2, 2.8284271247461900976L},
        {3, 14.646120160892684L},
        {5, 124.10555788745795L},
    };

    for (const auto& cp : checks) {
        const long double got = solve_for_p(cp.p, allow_multithreading, requested_threads, true);
        const long double err = std::fabsl(got - cp.expected);
        if (err > 5e-10L) {
            std::cerr << std::setprecision(18)
                      << "Checkpoint failed for P=" << cp.p << ": expected " << cp.expected
                      << ", got " << got << ", |err|=" << err << "\n";
            return false;
        }
    }

    std::cerr << "Validation checkpoints passed.\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints &&
        !run_checkpoints(options.allow_multithreading, options.requested_threads)) {
        return 1;
    }

    const long double ans =
        solve_for_p(options.p, options.allow_multithreading, options.requested_threads, false);

    std::cout << std::fixed << std::setprecision(10) << ans << '\n';
    return 0;
}
