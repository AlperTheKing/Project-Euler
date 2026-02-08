#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kDefaultLimit = 2000000000ULL;
constexpr u64 kMinGoodPrime = 193ULL;

struct Options {
    u64 limit = kDefaultLimit;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    bool run_checkpoints = true;
};

u64 isqrt_u64(const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / (x + 1ULL)) {
        ++x;
    }
    while (x > n / x) {
        --x;
    }
    return x;
}

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

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    long long parsed = 0;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10LL + static_cast<long long>(ch - '0');
        if (parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    value = static_cast<int>(parsed);
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
        if (parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.threads <= 0) {
        options.threads = 1;
    }

    return true;
}

std::vector<int> sieve_primes_up_to(const int n) {
    if (n < 2) {
        return {};
    }

    std::vector<char> is_prime(static_cast<std::size_t>(n) + 1U, 1);
    is_prime[0] = 0;
    is_prime[1] = 0;

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    for (int p = 2; p <= root; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        const int start = p * p;
        for (int x = start; x <= n; x += p) {
            is_prime[x] = 0;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / std::max(1.0, std::log(n))));
    for (int p = 2; p <= n; ++p) {
        if (is_prime[p]) {
            primes.push_back(p);
        }
    }
    return primes;
}

bool is_good_prime(const u64 p) {
    if (p <= 7ULL) {
        return false;
    }
    const u64 residue = p % 168ULL;
    return residue == 1ULL || residue == 25ULL || residue == 121ULL;
}

u64 singleton_sum_in_range(const u64 low,
                           const u64 high,
                           const u64 limit,
                           const std::vector<int>& base_primes) {
    if (low > high || high < 2ULL) {
        return 0ULL;
    }

    constexpr u64 segment_size = 1ULL << 20;
    u64 total = 0ULL;

    for (u64 segment_low = low; segment_low <= high;) {
        const u64 segment_high =
            std::min(segment_low + segment_size - 1ULL, high);
        const std::size_t length =
            static_cast<std::size_t>(segment_high - segment_low + 1ULL);

        std::vector<char> is_prime(length, 1);

        for (const int p_int : base_primes) {
            const u64 p = static_cast<u64>(p_int);
            const u64 p2 = p * p;
            if (p2 > segment_high) {
                break;
            }

            u64 start = (segment_low + p - 1ULL) / p;
            start *= p;
            if (start < p2) {
                start = p2;
            }

            for (u64 x = start; x <= segment_high; x += p) {
                is_prime[static_cast<std::size_t>(x - segment_low)] = 0;
            }
        }

        if (segment_low == 0ULL) {
            is_prime[0] = 0;
            if (length > 1U) {
                is_prime[1] = 0;
            }
        } else if (segment_low == 1ULL) {
            is_prime[0] = 0;
        }

        for (u64 value = std::max<u64>(2ULL, segment_low); value <= segment_high;
             ++value) {
            if (!is_prime[static_cast<std::size_t>(value - segment_low)]) {
                continue;
            }
            if (!is_good_prime(value)) {
                continue;
            }
            total += isqrt_u64(limit / value);
        }

        if (segment_high == std::numeric_limits<u64>::max()) {
            break;
        }
        segment_low = segment_high + 1ULL;
    }

    return total;
}

u64 sum_singleton_kernels(const u64 limit, const int threads) {
    if (limit < 11ULL) {
        return 0ULL;
    }

    const int base_bound = static_cast<int>(isqrt_u64(limit));
    const std::vector<int> base_primes = sieve_primes_up_to(base_bound);

    const u64 first = 2ULL;
    const u64 last = limit;
    const u64 count = last - first + 1ULL;

    int worker_count = std::max(1, threads);
    if (static_cast<u64>(worker_count) > count) {
        worker_count = static_cast<int>(count);
    }

    std::vector<u64> partial(static_cast<std::size_t>(worker_count), 0ULL);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(worker_count));

    const u64 block = (count + static_cast<u64>(worker_count) - 1ULL) /
                      static_cast<u64>(worker_count);

    for (int tid = 0; tid < worker_count; ++tid) {
        const u64 low = first + static_cast<u64>(tid) * block;
        if (low > last) {
            break;
        }
        const u64 high = std::min(last, low + block - 1ULL);

        workers.emplace_back([&, tid, low, high]() {
            partial[static_cast<std::size_t>(tid)] =
                singleton_sum_in_range(low, high, limit, base_primes);
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u64 total = 0ULL;
    for (const u64 value : partial) {
        total += value;
    }
    return total;
}

std::vector<int> collect_good_primes_up_to(const int n) {
    if (n < 11) {
        return {};
    }

    const std::vector<int> primes = sieve_primes_up_to(n);
    std::vector<int> good;
    good.reserve(primes.size() / 12U);

    for (const int p : primes) {
        if (is_good_prime(static_cast<u64>(p))) {
            good.push_back(p);
        }
    }

    return good;
}

u64 sum_multi_prime_kernels_dfs(const std::vector<int>& good_primes,
                                const u64 limit,
                                const std::size_t start_index,
                                const u64 current_product,
                                const int picked_count) {
    u64 total = 0ULL;

    for (std::size_t i = start_index; i < good_primes.size(); ++i) {
        const u64 p = static_cast<u64>(good_primes[i]);
        if (current_product > limit / p) {
            break;
        }

        const u64 next_product = current_product * p;
        const int next_picked = picked_count + 1;

        if (next_picked >= 2) {
            total += isqrt_u64(limit / next_product);
        }

        total += sum_multi_prime_kernels_dfs(good_primes,
                                             limit,
                                             i + 1,
                                             next_product,
                                             next_picked);
    }

    return total;
}

u64 sum_multi_prime_kernels(const u64 limit) {
    if (limit < kMinGoodPrime * kMinGoodPrime) {
        return 0ULL;
    }

    const u64 max_factor = limit / kMinGoodPrime;
    if (max_factor > static_cast<u64>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("max_factor does not fit into int");
    }

    const std::vector<int> good_primes =
        collect_good_primes_up_to(static_cast<int>(max_factor));

    return sum_multi_prime_kernels_dfs(good_primes, limit, 0U, 1ULL, 0);
}

u64 count_integer_representable(const u64 limit, const int threads) {
    const u64 square_kernel_contribution = isqrt_u64(limit);
    const u64 singleton_contribution = sum_singleton_kernels(limit, threads);
    const u64 multi_prime_contribution = sum_multi_prime_kernels(limit);

    return square_kernel_contribution + singleton_contribution +
           multi_prime_contribution;
}

std::vector<int> build_spf_table(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n) + 1U, 0);
    for (int i = 0; i <= n; ++i) {
        spf[static_cast<std::size_t>(i)] = i;
    }

    if (n >= 0) {
        spf[0] = 0;
    }
    if (n >= 1) {
        spf[1] = 1;
    }

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    for (int i = 2; i <= root; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (int x = i * i; x <= n; x += i) {
            if (spf[static_cast<std::size_t>(x)] == x) {
                spf[static_cast<std::size_t>(x)] = i;
            }
        }
    }

    return spf;
}

void add_factor(std::vector<std::pair<int, int>>& factors,
                const int prime,
                const int exponent) {
    for (auto& [p, e] : factors) {
        if (p == prime) {
            e += exponent;
            return;
        }
    }
    factors.push_back({prime, exponent});
}

std::vector<std::pair<int, int>> factorize_dy2(const int y,
                                                const int d,
                                                const std::vector<int>& spf) {
    std::vector<std::pair<int, int>> factors;

    int value = y;
    while (value > 1) {
        const int p = spf[static_cast<std::size_t>(value)];
        int exponent = 0;
        while (value % p == 0) {
            value /= p;
            ++exponent;
        }
        add_factor(factors, p, 2 * exponent);
    }

    int remaining_d = d;
    for (int p = 2; static_cast<long long>(p) * p <= remaining_d; ++p) {
        int exponent = 0;
        while (remaining_d % p == 0) {
            remaining_d /= p;
            ++exponent;
        }
        if (exponent > 0) {
            add_factor(factors, p, exponent);
        }
    }
    if (remaining_d > 1) {
        add_factor(factors, remaining_d, 1);
    }

    return factors;
}

void build_divisors_recursive(const std::vector<std::pair<int, int>>& factors,
                              const std::size_t index,
                              const u64 current,
                              std::vector<u64>& out) {
    if (index == factors.size()) {
        out.push_back(current);
        return;
    }

    const int prime = factors[index].first;
    const int exponent = factors[index].second;

    u64 value = current;
    for (int i = 0; i <= exponent; ++i) {
        build_divisors_recursive(factors, index + 1U, value, out);
        value *= static_cast<u64>(prime);
    }
}

std::vector<char> square_roots_with_positive_rep(const u64 limit,
                                                 const int d,
                                                 const std::vector<int>& spf) {
    const int max_root = static_cast<int>(isqrt_u64(limit));
    const int max_y = static_cast<int>(isqrt_u64(limit / static_cast<u64>(d)));

    std::vector<char> mark(static_cast<std::size_t>(max_root) + 1U, 0);

    std::vector<u64> divisors;
    divisors.reserve(512);

    for (int y = 1; y <= max_y; ++y) {
        const u64 D = static_cast<u64>(d) * static_cast<u64>(y) *
                      static_cast<u64>(y);

        const std::vector<std::pair<int, int>> factors = factorize_dy2(y, d, spf);

        divisors.clear();
        build_divisors_recursive(factors, 0U, 1ULL, divisors);

        for (const u64 u : divisors) {
            const u64 v = D / u;
            if (u > v) {
                continue;
            }
            if (((u + v) & 1ULL) != 0ULL) {
                continue;
            }
            if (v <= u) {
                continue;
            }

            const u64 m = (u + v) / 2ULL;
            if (m <= static_cast<u64>(max_root)) {
                mark[static_cast<std::size_t>(m)] = 1;
            }
        }
    }

    return mark;
}

u64 count_good_squares(const u64 limit) {
    const int max_root = static_cast<int>(isqrt_u64(limit));
    const std::vector<int> spf = build_spf_table(max_root);

    const std::vector<char> d1 = square_roots_with_positive_rep(limit, 1, spf);
    const std::vector<char> d2 = square_roots_with_positive_rep(limit, 2, spf);
    const std::vector<char> d3 = square_roots_with_positive_rep(limit, 3, spf);
    const std::vector<char> d7 = square_roots_with_positive_rep(limit, 7, spf);

    u64 count = 0ULL;
    for (int m = 1; m <= max_root; ++m) {
        if (d1[static_cast<std::size_t>(m)] && d2[static_cast<std::size_t>(m)] &&
            d3[static_cast<std::size_t>(m)] && d7[static_cast<std::size_t>(m)]) {
            ++count;
        }
    }
    return count;
}

u64 count_special_numbers(const u64 limit, const int threads) {
    const u64 integer_representable = count_integer_representable(limit, threads);
    const u64 total_squares = isqrt_u64(limit);
    const u64 good_squares = count_good_squares(limit);

    return integer_representable - total_squares + good_squares;
}

bool has_positive_rep_bruteforce(const u64 n, const int d) {
    for (u64 b = 1ULL;; ++b) {
        const u64 db2 = static_cast<u64>(d) * b * b;
        if (db2 >= n) {
            break;
        }

        const u64 rem = n - db2;
        const u64 a = isqrt_u64(rem);
        if (a > 0ULL && a * a == rem) {
            return true;
        }
    }
    return false;
}

u64 brute_force_count(const int limit) {
    const std::array<int, 4> ds = {1, 2, 3, 7};
    std::array<std::vector<char>, 4> mark;

    for (int i = 0; i < 4; ++i) {
        mark[i].assign(static_cast<std::size_t>(limit) + 1U, 0);
        const int d = ds[i];

        for (int b = 1;; ++b) {
            const u64 db2 = static_cast<u64>(d) * static_cast<u64>(b) *
                            static_cast<u64>(b);
            if (db2 > static_cast<u64>(limit)) {
                break;
            }
            const int max_a =
                static_cast<int>(isqrt_u64(static_cast<u64>(limit) - db2));
            for (int a = 1; a <= max_a; ++a) {
                const int n = static_cast<int>(static_cast<u64>(a) *
                                                   static_cast<u64>(a) +
                                               db2);
                mark[i][static_cast<std::size_t>(n)] = 1;
            }
        }
    }

    u64 count = 0ULL;
    for (int n = 1; n <= limit; ++n) {
        if (mark[0][static_cast<std::size_t>(n)] &&
            mark[1][static_cast<std::size_t>(n)] &&
            mark[2][static_cast<std::size_t>(n)] &&
            mark[3][static_cast<std::size_t>(n)]) {
            ++count;
        }
    }

    return count;
}

void run_checkpoints(const int threads) {
    const std::array<u64, 2> examples = {3600ULL, 88201ULL};
    for (const u64 n : examples) {
        for (const int d : {1, 2, 3, 7}) {
            if (!has_positive_rep_bruteforce(n, d)) {
                throw std::runtime_error("Example validation failed");
            }
        }
    }

    constexpr int small_limit = 200000;
    const u64 fast_small =
        count_special_numbers(static_cast<u64>(small_limit), threads);
    const u64 brute_small = brute_force_count(small_limit);
    if (fast_small != brute_small) {
        throw std::runtime_error("Small brute-force checkpoint failed");
    }

    constexpr u64 stated_limit = 10000000ULL;
    constexpr u64 stated_value = 75373ULL;
    const u64 fast_stated = count_special_numbers(stated_limit, threads);
    if (fast_stated != stated_value) {
        throw std::runtime_error("10^7 checkpoint failed");
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        std::cerr << "Usage: ./Euler229 [--limit=N] [--threads=T] [--skip-checkpoints]\n";
        return 1;
    }

    try {
        if (options.run_checkpoints) {
            run_checkpoints(options.threads);
        }

        const u64 answer = count_special_numbers(options.limit, options.threads);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
