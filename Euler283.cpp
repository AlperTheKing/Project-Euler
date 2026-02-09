#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i64 = std::int64_t;

constexpr int kDefaultLimit = 1000;

struct Options {
    int k_limit = kDefaultLimit;
    unsigned threads = std::thread::hardware_concurrency();
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
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg, const std::string& prefix, unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    unsigned parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10U + static_cast<unsigned>(c - '0');
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
        if (parse_int_after_prefix(arg, "--k-limit=", options.k_limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.k_limit <= 0) {
        return false;
    }
    if (options.threads == 0) {
        options.threads = 1;
    }
    return true;
}

std::vector<int> build_spf(const int limit) {
    std::vector<int> spf(static_cast<std::size_t>(limit + 1), 0);
    if (limit >= 1) {
        spf[1] = 1;
    }
    for (int i = 2; i <= limit; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        spf[static_cast<std::size_t>(i)] = i;
        if (static_cast<i64>(i) * i > limit) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            if (spf[static_cast<std::size_t>(j)] == 0) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return spf;
}

void factorize_int(int n, const std::vector<int>& spf, std::vector<std::pair<int, int>>& out) {
    out.clear();
    while (n > 1) {
        const int p = spf[static_cast<std::size_t>(n)];
        int cnt = 0;
        while (n % p == 0) {
            n /= p;
            ++cnt;
        }
        out.push_back({p, cnt});
    }
}

void merge_factorizations(const std::vector<std::pair<int, int>>& a,
                         const std::vector<std::pair<int, int>>& b,
                         std::vector<std::pair<int, int>>& out) {
    out.clear();
    std::size_t i = 0;
    std::size_t j = 0;

    while (i < a.size() || j < b.size()) {
        if (j == b.size() || (i < a.size() && a[i].first < b[j].first)) {
            out.push_back(a[i]);
            ++i;
        } else if (i == a.size() || b[j].first < a[i].first) {
            out.push_back(b[j]);
            ++j;
        } else {
            out.push_back({a[i].first, a[i].second + b[j].second});
            ++i;
            ++j;
        }
    }
}

void generate_divisors(const std::vector<std::pair<int, int>>& factors, std::vector<u64>& divisors) {
    divisors.clear();
    divisors.push_back(1);

    for (const auto& [p, e] : factors) {
        const std::size_t base = divisors.size();
        u64 mul = 1;
        for (int i = 1; i <= e; ++i) {
            mul *= static_cast<u64>(p);
            for (std::size_t j = 0; j < base; ++j) {
                divisors.push_back(divisors[j] * mul);
            }
        }
    }
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 solve_range(const int k_start, const int k_end, const std::vector<int>& spf) {
    std::vector<std::pair<int, int>> factors_n;
    std::vector<std::pair<int, int>> factors_m;
    std::vector<std::pair<int, int>> factors_all;
    std::vector<u64> divisors;

    u128 total = 0;

    for (int k = k_start; k <= k_end; ++k) {
        const u64 r = static_cast<u64>(2 * k);
        const u64 n = r * r;
        factorize_int(static_cast<int>(n), spf, factors_n);

        const u64 x_max = static_cast<u64>(std::sqrt(static_cast<long double>(3 * n)));
        for (u64 x = 1; x <= x_max; ++x) {
            const u64 m = x * x + n;
            factorize_int(static_cast<int>(m), spf, factors_m);
            merge_factorizations(factors_n, factors_m, factors_all);
            generate_divisors(factors_all, divisors);

            const u64 M = n * m;
            for (u64 u : divisors) {
                if (u * u > M) {
                    continue;
                }

                const u64 v = M / u;
                if ((u + n) % x != 0 || (v + n) % x != 0) {
                    continue;
                }

                const u64 y = (u + n) / x;
                if (y < x) {
                    continue;
                }
                const u64 z = (v + n) / x;

                total += static_cast<u128>(2) * (x + y + z);
            }
        }
    }

    return total;
}

u128 solve_fast(const int k_limit, const std::vector<int>& spf, unsigned threads) {
    if (threads <= 1 || k_limit < 100) {
        return solve_range(1, k_limit, spf);
    }

    const unsigned use_threads = std::min<unsigned>(threads, static_cast<unsigned>(k_limit));
    const int chunk = (k_limit + static_cast<int>(use_threads) - 1) / static_cast<int>(use_threads);

    std::vector<std::thread> pool;
    std::vector<u128> partial(use_threads, 0);
    pool.reserve(use_threads);

    for (unsigned t = 0; t < use_threads; ++t) {
        const int start = static_cast<int>(t) * chunk + 1;
        const int end = std::min(k_limit, start + chunk - 1);
        if (start > end) {
            continue;
        }

        pool.emplace_back([&, start, end, t]() {
            partial[t] = solve_range(start, end, spf);
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u128 total = 0;
    for (u128 part : partial) {
        total += part;
    }
    return total;
}

u128 solve_bruteforce_small(const int k_limit) {
    u128 total = 0;

    for (int k = 1; k <= k_limit; ++k) {
        const u64 r = static_cast<u64>(2 * k);
        const u64 n = r * r;

        const u64 x_max = static_cast<u64>(std::sqrt(static_cast<long double>(3 * n)));
        for (u64 x = 1; x <= x_max; ++x) {
            const long double disc = std::sqrt(static_cast<long double>(n) * (static_cast<long double>(n) + x * x));
            const u64 y_max = static_cast<u64>((n + disc) / x + 1e-12L);

            for (u64 y = x; y <= y_max; ++y) {
                const i64 d = static_cast<i64>(x * y) - static_cast<i64>(n);
                if (d <= 0) {
                    continue;
                }

                const u64 num = n * (x + y);
                if (num % static_cast<u64>(d) != 0) {
                    continue;
                }

                const u64 z = num / static_cast<u64>(d);
                if (z < y) {
                    continue;
                }

                if (x * y * z != n * (x + y + z)) {
                    continue;
                }

                total += static_cast<u128>(2) * (x + y + z);
            }
        }
    }

    return total;
}

bool run_checkpoints(const std::vector<int>& spf) {
    const u128 fast5 = solve_fast(5, spf, 1);
    const u128 brute5 = solve_bruteforce_small(5);
    if (fast5 != brute5 || fast5 != static_cast<u128>(140098ULL)) {
        std::cerr << "Checkpoint failed for k<=5: fast=" << to_string_u128(fast5)
                  << ", brute=" << to_string_u128(brute5) << '\n';
        return false;
    }

    const u128 fast10 = solve_fast(10, spf, 1);
    const u128 brute10 = solve_bruteforce_small(10);
    if (fast10 != brute10 || fast10 != static_cast<u128>(3781786ULL)) {
        std::cerr << "Checkpoint failed for k<=10: fast=" << to_string_u128(fast10)
                  << ", brute=" << to_string_u128(brute10) << '\n';
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

    // max(x^2 + n) with n=(2k)^2 and x<=sqrt(3n) is at most 4*n
    const int n_max = 4 * (2 * options.k_limit) * (2 * options.k_limit);
    const std::vector<int> spf = build_spf(n_max);

    if (options.run_checkpoints && !run_checkpoints(spf)) {
        return 2;
    }

    const u128 answer = solve_fast(options.k_limit, spf, options.threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
