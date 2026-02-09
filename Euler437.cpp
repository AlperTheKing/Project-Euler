#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 100000000;
    int threads = 0;
    bool run_checkpoints = true;
};

struct Result {
    u64 sum = 0;
    int count = 0;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2 && options.threads >= 0;
}

int choose_thread_count(int requested, std::size_t work_items) {
    if (work_items <= 1) {
        return 1;
    }
    int threads = requested;
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
    }
    if (threads <= 0) {
        threads = 4;
    }
    if (threads > static_cast<int>(work_items)) {
        threads = static_cast<int>(work_items);
    }
    if (threads < 1) {
        threads = 1;
    }
    return threads;
}

std::vector<std::uint16_t> build_odd_spf(const int limit) {
    std::vector<std::uint16_t> spf(static_cast<std::size_t>(limit / 2 + 1), 0);
    const int root = static_cast<int>(std::sqrt(static_cast<double>(limit)));
    for (int i = 3; i <= root; i += 2) {
        if (spf[static_cast<std::size_t>(i >> 1)] != 0) {
            continue;
        }
        const int step = i << 1;
        for (int j = i * i; j <= limit; j += step) {
            std::uint16_t& slot = spf[static_cast<std::size_t>(j >> 1)];
            if (slot == 0) {
                slot = static_cast<std::uint16_t>(i);
            }
        }
    }
    return spf;
}

std::vector<int> collect_primes(const int limit, const std::vector<std::uint16_t>& spf) {
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / 12));
    if (limit >= 2) {
        primes.push_back(2);
    }
    for (int p = 3; p <= limit; p += 2) {
        if (spf[static_cast<std::size_t>(p >> 1)] == 0) {
            primes.push_back(p);
        }
    }
    return primes;
}

inline u64 mul_mod(u64 a, u64 b, int mod) {
    return (a * b) % static_cast<u64>(mod);
}

int mod_pow(int base, u64 exp, int mod) {
    u64 result = 1ULL;
    u64 cur = static_cast<u64>(base % mod + mod) % static_cast<u64>(mod);
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = mul_mod(result, cur, mod);
        }
        cur = mul_mod(cur, cur, mod);
        e >>= 1ULL;
    }
    return static_cast<int>(result);
}

int tonelli_sqrt_5(int p) {
    if (p == 2) {
        return 1;
    }
    if (mod_pow(5, static_cast<u64>((p - 1) / 2), p) != 1) {
        return -1;
    }
    if ((p & 3) == 3) {
        return mod_pow(5, static_cast<u64>((p + 1) / 4), p);
    }

    int q = p - 1;
    int s = 0;
    while ((q & 1) == 0) {
        q >>= 1;
        ++s;
    }

    int z = 2;
    while (mod_pow(z, static_cast<u64>((p - 1) / 2), p) != p - 1) {
        ++z;
    }

    u64 c = static_cast<u64>(mod_pow(z, static_cast<u64>(q), p));
    u64 x = static_cast<u64>(mod_pow(5, static_cast<u64>((q + 1) / 2), p));
    u64 t = static_cast<u64>(mod_pow(5, static_cast<u64>(q), p));
    int m = s;

    while (t != 1ULL) {
        int i = 1;
        u64 tt = mul_mod(t, t, p);
        while (i < m && tt != 1ULL) {
            tt = mul_mod(tt, tt, p);
            ++i;
        }

        const int shift = m - i - 1;
        const u64 b = static_cast<u64>(mod_pow(static_cast<int>(c), 1ULL << shift, p));
        x = mul_mod(x, b, p);
        t = mul_mod(t, mul_mod(b, b, p), p);
        c = mul_mod(b, b, p);
        m = i;
    }

    return static_cast<int>(x);
}

int factor_unique(int n, const std::vector<std::uint16_t>& spf, std::array<int, 12>& out) {
    int cnt = 0;
    int x = n;

    if ((x & 1) == 0) {
        out[static_cast<std::size_t>(cnt++)] = 2;
        while ((x & 1) == 0) {
            x >>= 1;
        }
    }

    while (x > 1) {
        int f = 0;
        if ((x & 1) == 1) {
            f = static_cast<int>(spf[static_cast<std::size_t>(x >> 1)]);
            if (f == 0) {
                f = x;
            }
        } else {
            f = 2;
        }

        out[static_cast<std::size_t>(cnt++)] = f;
        while (x % f == 0) {
            x /= f;
        }
    }

    return cnt;
}

bool is_primitive_root(int g, int p, const std::array<int, 12>& factors, int factor_count) {
    const int phi = p - 1;
    for (int i = 0; i < factor_count; ++i) {
        const int q = factors[static_cast<std::size_t>(i)];
        if (mod_pow(g, static_cast<u64>(phi / q), p) == 1) {
            return false;
        }
    }
    return true;
}

Result solve(int limit, int requested_threads) {
    const std::vector<std::uint16_t> spf = build_odd_spf(limit);
    const std::vector<int> primes = collect_primes(limit, spf);

    const int thread_count = choose_thread_count(requested_threads, primes.size());
    std::vector<u64> partial_sum(static_cast<std::size_t>(thread_count), 0ULL);
    std::vector<int> partial_count(static_cast<std::size_t>(thread_count), 0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        const std::size_t begin =
            primes.size() * static_cast<std::size_t>(t) / static_cast<std::size_t>(thread_count);
        const std::size_t end =
            primes.size() * static_cast<std::size_t>(t + 1) / static_cast<std::size_t>(thread_count);

        workers.emplace_back([&, begin, end, t]() {
            u64 local_sum = 0ULL;
            int local_count = 0;
            std::array<int, 12> factors{};

            for (std::size_t idx = begin; idx < end; ++idx) {
                const int p = primes[idx];
                if (p == 2 || p == 3) {
                    continue;
                }
                if (p == 5) {
                    local_sum += 5ULL;
                    ++local_count;
                    continue;
                }

                const int mod5 = p % 5;
                if (mod5 != 1 && mod5 != 4) {
                    continue;
                }

                const int sqrt5 = tonelli_sqrt_5(p);
                if (sqrt5 < 0) {
                    continue;
                }

                const int inv2 = (p + 1) / 2;
                const int r1 = static_cast<int>((static_cast<u64>((1 + sqrt5) % p) * inv2) % p);
                const int r2 = static_cast<int>((static_cast<u64>((1 + p - sqrt5) % p) * inv2) % p);

                const int fc = factor_unique(p - 1, spf, factors);
                if (is_primitive_root(r1, p, factors, fc) || is_primitive_root(r2, p, factors, fc)) {
                    local_sum += static_cast<u64>(p);
                    ++local_count;
                }
            }

            partial_sum[static_cast<std::size_t>(t)] = local_sum;
            partial_count[static_cast<std::size_t>(t)] = local_count;
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    Result res;
    for (int t = 0; t < thread_count; ++t) {
        res.sum += partial_sum[static_cast<std::size_t>(t)];
        res.count += partial_count[static_cast<std::size_t>(t)];
    }
    return res;
}

bool run_checkpoints(int requested_threads) {
    const Result small = solve(10000, requested_threads);
    if (small.count != 323) {
        std::cerr << "Checkpoint failed: count below 10000\n";
        return false;
    }
    if (small.sum != 1480491ULL) {
        std::cerr << "Checkpoint failed: sum below 10000\n";
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
    if (options.run_checkpoints && !run_checkpoints(options.threads)) {
        return 2;
    }

    const Result ans = solve(options.limit, options.threads);
    std::cout << ans.sum << '\n';
    return 0;
}
