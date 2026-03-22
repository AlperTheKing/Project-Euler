#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <string>
#include <thread>
#include <vector>

namespace {

using i128 = __int128_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 MOD = 1'000'000'009ULL;
constexpr u64 TARGET_LIMIT = 100'000'000'000'000ULL;
constexpr u64 SAMPLE_LIMIT = 1'000ULL;
constexpr u64 SAMPLE_SUM = 190'950'976ULL;

constexpr u64 SQRT5_MOD = 383'008'016ULL;
constexpr u64 PHI_MOD = 691'504'013ULL;
constexpr u64 PSI_MOD = 308'495'997ULL;

struct Options {
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
};

void usage() {
    std::cerr
        << "Usage:\n"
        << "  ./Euler989 validate [check_max] [--single-thread] [--threads=N]\n"
        << "  ./Euler989 sum <limit> [--single-thread] [--threads=N]\n"
        << "  ./Euler989 answer [--single-thread] [--threads=N]\n";
}

u64 mod_mul(const u64 a, const u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(MOD));
}

u64 mod_add(const u64 a, const u64 b) {
    const u64 s = a + b;
    return s >= MOD ? s - MOD : s;
}

u64 mod_sub(const u64 a, const u64 b) {
    return a >= b ? a - b : a + MOD - b;
}

u64 mod_neg(const u64 a) {
    return a == 0 ? 0 : MOD - a;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        exp >>= 1U;
    }
    return result;
}

u64 isqrt(const u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }
    while (x > n / x) {
        --x;
    }
    return x;
}

std::vector<u32> sieve_primes(const int limit) {
    if (limit < 2) {
        return {};
    }
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    std::vector<u32> primes;
    is_prime[0] = false;
    is_prime[1] = false;
    for (int i = 2; i <= limit; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        primes.push_back(static_cast<u32>(i));
        if (i > limit / i) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }
    return primes;
}

std::vector<std::int8_t> mobius_sieve(const int limit) {
    std::vector<std::int8_t> mu(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> primes;
    std::vector<int> least(static_cast<std::size_t>(limit + 1), 0);
    mu[1] = 1;

    for (int i = 2; i <= limit; ++i) {
        if (least[static_cast<std::size_t>(i)] == 0) {
            least[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (const int p : primes) {
            const int ip = i * p;
            if (ip > limit || p > least[static_cast<std::size_t>(i)]) {
                break;
            }
            least[static_cast<std::size_t>(ip)] = p;
            if (p == least[static_cast<std::size_t>(i)]) {
                mu[static_cast<std::size_t>(ip)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(ip)] = -mu[static_cast<std::size_t>(i)];
        }
    }

    return mu;
}

u32 g_bruteforce(const u64 n) {
    u32 count = 0;
    for (u64 x = 0; x < n; ++x) {
        if ((x * x + n - x - 1) % n == 0) {
            ++count;
        }
    }
    return count;
}

u32 g_from_factorization(u64 n, const std::vector<u32>& primes) {
    if (n == 1) {
        return 1;
    }

    u32 result = 1;
    for (const u32 p32 : primes) {
        const u64 p = static_cast<u64>(p32);
        if (p > n / p) {
            break;
        }
        if (n % p != 0) {
            continue;
        }

        u32 exponent = 0;
        while (n % p == 0) {
            n /= p;
            ++exponent;
        }

        if (p == 2) {
            return 0;
        }
        if (p == 5) {
            if (exponent >= 2) {
                return 0;
            }
            continue;
        }

        switch (p % 5) {
            case 1:
            case 4:
                result = static_cast<u32>(result * 2U);
                break;
            case 2:
            case 3:
                return 0;
            default:
                std::abort();
        }
    }

    if (n > 1) {
        if (n == 2) {
            return 0;
        }
        if (n == 5) {
            return result;
        }
        switch (n % 5) {
            case 1:
            case 4:
                result = static_cast<u32>(result * 2U);
                break;
            case 2:
            case 3:
                return 0;
            default:
                std::abort();
        }
    }

    return result;
}

u64 q_form(const u64 a, const u64 b) {
    return a * a - a * b - b * b;
}

u64 gcd(u64 a, u64 b) {
    while (b != 0) {
        const u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u32 reduced_pair_count(const u64 n) {
    u32 count = 0;
    const u64 max_b = isqrt(n);
    for (u64 b = 1; b <= max_b; ++b) {
        for (u64 a = 2 * b;; ++a) {
            const u64 q = q_form(a, b);
            if (q > n) {
                break;
            }
            if (q == n && gcd(a, b) == 1) {
                ++count;
            }
        }
    }
    return count;
}

u64 direct_pair_sum(const u64 limit, const u64 base) {
    u64 acc = 0;
    const u64 max_b = isqrt(limit);
    for (u64 b = 1; b <= max_b; ++b) {
        for (u64 a = 2 * b;; ++a) {
            const u64 q = q_form(a, b);
            if (q > limit) {
                break;
            }
            if (gcd(a, b) == 1) {
                acc = mod_add(acc, mod_pow(base, q));
            }
        }
    }
    return acc;
}

struct PowerSeq {
    u64 v = 0;
    u64 term = 1;
    u64 delta = 1;
    u64 delta_step = 1;

    static PowerSeq square(const u64 base) {
        const u64 base_sq = mod_mul(base, base);
        return {0, 1, base, base_sq};
    }

    static PowerSeq triangular(const u64 base) {
        const u64 base_sq = mod_mul(base, base);
        return {0, 1, base_sq, base_sq};
    }

    void step() {
        term = mod_mul(term, delta);
        delta = mod_mul(delta, delta_step);
        ++v;
    }

    void extend_through(const u64 target, u64& window) {
        while (v <= target) {
            window = mod_add(window, term);
            step();
        }
    }

    void trim_before(const u64 target, u64& window) {
        while (v < target) {
            window = mod_sub(window, term);
            step();
        }
    }
};

u64 nonprimitive_sum(const u64 limit, const u64 w, const u64 winv) {
    if (limit == 0) {
        return 0;
    }

    const u64 sqrt_limit = isqrt(limit);
    const u64 winv5 = mod_pow(winv, 5);
    const u64 winv10 = mod_mul(winv5, winv5);
    u64 ans = 0;

    const u64 even_t_max = sqrt_limit / 2;
    if (even_t_max > 0) {
        PowerSeq add_seq = PowerSeq::square(w);
        PowerSeq trim_seq = PowerSeq::square(w);
        u64 window = 0;
        u64 factor = 1;
        u64 ratio = winv5;

        for (u64 t = 1; t <= even_t_max; ++t) {
            factor = mod_mul(factor, ratio);
            ratio = mod_mul(ratio, winv10);

            const u64 vmax = isqrt(limit + 5 * t * t);
            add_seq.extend_through(vmax, window);
            trim_seq.trim_before(3 * t, window);
            ans = mod_add(ans, mod_mul(factor, window));
        }
    }

    const u64 odd_t_max = (sqrt_limit - 1) / 2;
    PowerSeq add_seq = PowerSeq::triangular(w);
    PowerSeq trim_seq = PowerSeq::triangular(w);
    u64 window = 0;
    u64 factor = winv;
    u64 ratio = winv10;

    for (u64 t = 0; t <= odd_t_max; ++t) {
        const u64 disc = 4 * limit + 20 * t * t + 20 * t + 5;
        const u64 vmax = (isqrt(disc) - 1) / 2;
        add_seq.extend_through(vmax, window);
        trim_seq.trim_before(3 * t + 1, window);
        ans = mod_add(ans, mod_mul(factor, window));
        factor = mod_mul(factor, ratio);
        ratio = mod_mul(ratio, winv10);
    }

    return ans;
}

std::vector<u32> square_powers(const u64 base, const int max_g) {
    std::vector<u32> out(static_cast<std::size_t>(max_g + 1), 0);
    out[0] = 1;
    if (max_g == 0) {
        return out;
    }

    const u64 base_sq = mod_mul(base, base);
    u64 value = 1;
    u64 ratio = base;
    for (int g = 1; g <= max_g; ++g) {
        value = mod_mul(value, ratio);
        out[static_cast<std::size_t>(g)] = static_cast<u32>(value);
        ratio = mod_mul(ratio, base_sq);
    }
    return out;
}

u64 normalize_signed_mod(i128 value) {
    value %= static_cast<i128>(MOD);
    if (value < 0) {
        value += static_cast<i128>(MOD);
    }
    return static_cast<u64>(value);
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const int workload) {
    if (!allow_multithreading || workload <= 1) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }
    if (threads > static_cast<unsigned>(workload)) {
        threads = static_cast<unsigned>(workload);
    }
    return std::max(1U, threads);
}

struct PrimitiveWorkerTask {
    u64 limit = 0;
    int start_g = 1;
    int end_g = 1;
    const std::vector<std::int8_t>* mu = nullptr;
    const std::vector<u32>* phi_sq = nullptr;
    const std::vector<u32>* phi_inv_sq = nullptr;
    i128 phi_total = 0;
    i128 psi_total = 0;
};

void* primitive_worker_entry(void* arg) {
    auto& task = *static_cast<PrimitiveWorkerTask*>(arg);
    for (int g = task.start_g; g < task.end_g; ++g) {
        const int mu_g = static_cast<int>((*task.mu)[static_cast<std::size_t>(g)]);
        if (mu_g == 0) {
            continue;
        }

        const u64 gg = static_cast<u64>(g) * static_cast<u64>(g);
        const u64 scaled_limit = task.limit / gg;

        const u64 phi_w = static_cast<u64>((*task.phi_sq)[static_cast<std::size_t>(g)]);
        const u64 phi_winv = static_cast<u64>((*task.phi_inv_sq)[static_cast<std::size_t>(g)]);
        const u64 psi_w = (g & 1) == 0 ? phi_winv : mod_neg(phi_winv);
        const u64 psi_winv = (g & 1) == 0 ? phi_w : mod_neg(phi_w);

        const i128 phi_value = static_cast<i128>(nonprimitive_sum(scaled_limit, phi_w, phi_winv));
        const i128 psi_value = static_cast<i128>(nonprimitive_sum(scaled_limit, psi_w, psi_winv));

        if (mu_g > 0) {
            task.phi_total += phi_value;
            task.psi_total += psi_value;
        } else {
            task.phi_total -= phi_value;
            task.psi_total -= psi_value;
        }
    }
    return nullptr;
}

std::pair<u64, u64> primitive_sums(const u64 limit,
                                   const std::vector<std::int8_t>& mu,
                                   const std::vector<u32>& phi_sq,
                                   const std::vector<u32>& phi_inv_sq,
                                   const Options& options) {
    const int max_g = static_cast<int>(isqrt(limit));
    const unsigned thread_count = choose_thread_count(options.allow_multithreading,
                                                      options.requested_threads,
                                                      max_g);

    std::vector<pthread_t> threads(static_cast<std::size_t>(thread_count));
    std::vector<PrimitiveWorkerTask> tasks(static_cast<std::size_t>(thread_count));

    for (unsigned t = 0; t < thread_count; ++t) {
        PrimitiveWorkerTask& task = tasks[static_cast<std::size_t>(t)];
        task.limit = limit;
        task.start_g = 1 + static_cast<int>((static_cast<u64>(max_g) * t) / thread_count);
        task.end_g = 1 + static_cast<int>((static_cast<u64>(max_g) * (t + 1U)) / thread_count);
        task.mu = &mu;
        task.phi_sq = &phi_sq;
        task.phi_inv_sq = &phi_inv_sq;
        task.phi_total = 0;
        task.psi_total = 0;
    }

    for (unsigned t = 0; t < thread_count; ++t) {
        const int rc = pthread_create(&threads[static_cast<std::size_t>(t)],
                                      nullptr,
                                      primitive_worker_entry,
                                      &tasks[static_cast<std::size_t>(t)]);
        assert(rc == 0);
    }
    for (unsigned t = 0; t < thread_count; ++t) {
        const int rc = pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
    }

    i128 phi_total = 0;
    i128 psi_total = 0;
    for (const PrimitiveWorkerTask& task : tasks) {
        phi_total += task.phi_total;
        psi_total += task.psi_total;
    }

    return {normalize_signed_mod(phi_total), normalize_signed_mod(psi_total)};
}

u64 solve(const u64 limit, const Options& options) {
    const int max_g = static_cast<int>(isqrt(limit));
    const std::vector<std::int8_t> mu = mobius_sieve(max_g);

    const u64 phi_inv = mod_pow(PHI_MOD, MOD - 2);
    const std::vector<u32> phi_sq = square_powers(PHI_MOD, max_g);
    const std::vector<u32> phi_inv_sq = square_powers(phi_inv, max_g);

    const auto [phi_sum, psi_sum] = primitive_sums(limit, mu, phi_sq, phi_inv_sq, options);
    const u64 inv_sqrt5 = mod_pow(SQRT5_MOD, MOD - 2);
    return mod_mul(mod_sub(phi_sum, psi_sum), inv_sqrt5);
}

u64 checksum_via_factorization(const u64 limit, const std::vector<u32>& primes) {
    u64 acc = 0;
    u64 f_prev = 0;
    u64 f_cur = 1;
    for (u64 n = 1; n <= limit; ++n) {
        const u64 g = static_cast<u64>(g_from_factorization(n, primes));
        acc = mod_add(acc, mod_mul(f_cur, g));
        const u64 f_next = mod_add(f_prev, f_cur);
        f_prev = f_cur;
        f_cur = f_next;
    }
    return acc;
}

void validate(const u64 check_max, const Options& options) {
    assert(mod_mul(SQRT5_MOD, SQRT5_MOD) == 5);
    assert(mod_sub(mod_mul(PHI_MOD, PHI_MOD), PHI_MOD) == 1);
    assert(mod_sub(mod_mul(PSI_MOD, PSI_MOD), PSI_MOD) == 1);
    assert(mod_mul(PHI_MOD, PSI_MOD) == MOD - 1);
    std::cout << "Checkpoint 1 passed: modular sqrt(5), phi, and psi are consistent.\n";

    const int prime_limit = static_cast<int>(isqrt(check_max)) + 10;
    const std::vector<u32> primes = sieve_primes(prime_limit);
    for (u64 n = 1; n <= check_max; ++n) {
        const u32 brute = g_bruteforce(n);
        const u32 factorized = g_from_factorization(n, primes);
        assert(brute == factorized);
    }
    std::cout << "Checkpoint 2 passed: G(n) factorization matches brute force.\n";

    for (u64 n = 1; n <= check_max; ++n) {
        const u32 reduced_pairs = reduced_pair_count(n);
        const u32 factorized = g_from_factorization(n, primes);
        assert(reduced_pairs == factorized);
    }
    std::cout << "Checkpoint 3 passed: reduced quadratic-form pairs match G(n).\n";

    const u64 direct_limit = std::min<u64>(check_max, 250);
    const u64 phi_pair_sum = direct_pair_sum(direct_limit, PHI_MOD);
    const u64 psi_pair_sum = direct_pair_sum(direct_limit, PSI_MOD);
    const int max_g = static_cast<int>(isqrt(direct_limit));
    const std::vector<std::int8_t> mu = mobius_sieve(max_g);
    const u64 phi_inv = mod_pow(PHI_MOD, MOD - 2);
    const std::vector<u32> phi_sq = square_powers(PHI_MOD, max_g);
    const std::vector<u32> phi_inv_sq = square_powers(phi_inv, max_g);
    const auto [phi_fast, psi_fast] = primitive_sums(direct_limit, mu, phi_sq, phi_inv_sq, options);
    assert(phi_fast == phi_pair_sum);
    assert(psi_fast == psi_pair_sum);
    std::cout << "Checkpoint 4 passed: Möbius/parity solver matches direct pair enumeration.\n";

    const int sample_prime_limit = static_cast<int>(isqrt(SAMPLE_LIMIT)) + 10;
    const std::vector<u32> sample_primes = sieve_primes(sample_prime_limit);
    const u64 sample_fast = solve(SAMPLE_LIMIT, options);
    const u64 sample_factorized = checksum_via_factorization(SAMPLE_LIMIT, sample_primes);
    assert(sample_fast == SAMPLE_SUM);
    assert(sample_fast == sample_factorized);
    std::cout << "Checkpoint 5 passed: sample checksum equals " << sample_fast << ".\n";
}

bool parse_u64(const std::string& text, u64& value) {
    if (text.empty()) {
        return false;
    }
    u64 parsed = 0;
    for (const char c : text) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg, const char* prefix, unsigned& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }
    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    if (!parse_u64(tail, parsed) || parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_command_options(int argc,
                           char** argv,
                           int start_index,
                           Options& options,
                           std::vector<std::string>& positional) {
    for (int i = start_index; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        unsigned threads = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", threads)) {
            options.requested_threads = threads;
            continue;
        }
        positional.push_back(arg);
    }
    return true;
}

double elapsed_seconds(const std::chrono::steady_clock::time_point started) {
    const auto elapsed = std::chrono::steady_clock::now() - started;
    return std::chrono::duration<double>(elapsed).count();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        usage();
        return 0;
    }

    const std::string command(argv[1]);
    Options options;
    std::vector<std::string> positional;
    if (!parse_command_options(argc, argv, 2, options, positional)) {
        usage();
        return 1;
    }

    if (command == "validate") {
        if (positional.size() > 1) {
            usage();
            return 1;
        }
        u64 check_max = 200;
        if (!positional.empty() && !parse_u64(positional[0], check_max)) {
            usage();
            return 1;
        }
        const auto started = std::chrono::steady_clock::now();
        validate(check_max, options);
        std::cout << std::fixed << std::setprecision(3)
                  << "Validation completed in " << elapsed_seconds(started) << "s.\n";
        return 0;
    }

    if (command == "sum") {
        if (positional.size() != 1) {
            usage();
            return 1;
        }
        u64 limit = 0;
        if (!parse_u64(positional[0], limit)) {
            usage();
            return 1;
        }
        const auto started = std::chrono::steady_clock::now();
        const u64 answer = solve(limit, options);
        std::cout << answer << '\n';
        std::cerr << std::fixed << std::setprecision(3)
                  << "Computed in " << elapsed_seconds(started) << "s.\n";
        return 0;
    }

    if (command == "answer") {
        if (!positional.empty()) {
            usage();
            return 1;
        }
        const auto started = std::chrono::steady_clock::now();
        const u64 answer = solve(TARGET_LIMIT, options);
        std::cout << answer << '\n';
        std::cerr << std::fixed << std::setprecision(3)
                  << "Computed in " << elapsed_seconds(started) << "s.\n";
        return 0;
    }

    usage();
    return 0;
}
