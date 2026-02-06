#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u16 = std::uint16_t;
using u8 = std::uint8_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultLimit = 10'000'000'000ULL;
constexpr u64 kCheckpointKnownLimit = 100ULL;
constexpr u64 kCheckpointKnownExpected = 193ULL;
constexpr u64 kCheckpointBruteLimit = 500ULL;
constexpr u64 kThreadConsistencyLimit = 200'000ULL;

struct Options {
    u64 limit = kDefaultLimit;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
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

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--limit=", parsed_u64)) {
            options.limit = parsed_u64;
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

    if (options.limit < 2ULL) {
        std::cerr << "--limit must be at least 2.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2ULL) {
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

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 pow_u64_capped(u64 base, int exp, u64 cap) {
    u64 result = 1ULL;
    for (int i = 0; i < exp; ++i) {
        if (base != 0ULL && result > cap / base) {
            return cap + 1ULL;
        }
        result *= base;
    }
    return result;
}

u64 kth_root_floor(u64 n, int k) {
    if (k <= 1 || n <= 1ULL) {
        return n;
    }

    const long double approx =
        std::pow(static_cast<long double>(n), 1.0L / static_cast<long double>(k));
    u64 r = static_cast<u64>(approx);
    if (r == 0ULL) {
        r = 1ULL;
    }

    while (pow_u64_capped(r + 1ULL, k, n) <= n) {
        ++r;
    }
    while (pow_u64_capped(r, k, n) > n) {
        --r;
    }
    return r;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string s;
    while (value > 0) {
        const u32 digit = static_cast<u32>(value % 10U);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct PrimeTable {
    int limit = 0;
    std::vector<int> primes;
    std::vector<int> pi_small;

    explicit PrimeTable(int sieve_limit) { build(sieve_limit); }

    void build(int sieve_limit) {
        limit = sieve_limit;
        std::vector<int> lp(static_cast<std::size_t>(limit) + 1ULL, 0);
        primes.reserve(static_cast<std::size_t>(limit) / 10ULL + 32ULL);

        for (int i = 2; i <= limit; ++i) {
            if (lp[static_cast<std::size_t>(i)] == 0) {
                lp[static_cast<std::size_t>(i)] = i;
                primes.push_back(i);
            }

            for (const int p : primes) {
                const long long v = 1LL * p * i;
                if (p > lp[static_cast<std::size_t>(i)] || v > limit) {
                    break;
                }
                lp[static_cast<std::size_t>(v)] = p;
            }
        }

        pi_small.assign(static_cast<std::size_t>(limit) + 1ULL, 0);
        int ptr = 0;
        for (int i = 1; i <= limit; ++i) {
            pi_small[static_cast<std::size_t>(i)] =
                pi_small[static_cast<std::size_t>(i - 1)];
            if (ptr < static_cast<int>(primes.size()) && primes[static_cast<std::size_t>(ptr)] == i) {
                ++pi_small[static_cast<std::size_t>(i)];
                ++ptr;
            }
        }
    }
};

u64 icbrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while (static_cast<u128>(r + 1ULL) * (r + 1ULL) * (r + 1ULL) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * r * r > x) {
        --r;
    }
    return r;
}

u64 iroot4_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(std::sqrt(static_cast<long double>(x))));
    while (static_cast<u128>(r + 1ULL) * (r + 1ULL) * (r + 1ULL) * (r + 1ULL) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * r * r * r > x) {
        --r;
    }
    return r;
}

class PrimeCounter {
public:
    explicit PrimeCounter(const PrimeTable& table) : table_(table) { build_phi_table(); }

    u64 pi(u64 x) const {
        if (x <= static_cast<u64>(table_.limit)) {
            return static_cast<u64>(table_.pi_small[static_cast<std::size_t>(x)]);
        }

        auto& cache = pi_cache();
        const auto it = cache.find(x);
        if (it != cache.end()) {
            return it->second;
        }

        const u64 a = pi(iroot4_u64(x));
        const u64 b = pi(isqrt_u64(x));
        const u64 c = pi(icbrt_u64(x));

        i64 sum = static_cast<i64>(phi(x, static_cast<int>(a)));
        sum += static_cast<i64>((b + a - 2ULL) * (b - a + 1ULL) / 2ULL);

        for (u64 i = a + 1ULL; i <= b; ++i) {
            const u64 w = x / static_cast<u64>(table_.primes[static_cast<std::size_t>(i - 1ULL)]);
            sum -= static_cast<i64>(pi(w));
            if (i <= c) {
                const u64 bi = pi(isqrt_u64(w));
                for (u64 j = i; j <= bi; ++j) {
                    const u64 q =
                        w / static_cast<u64>(table_.primes[static_cast<std::size_t>(j - 1ULL)]);
                    sum -= static_cast<i64>(pi(q)) - static_cast<i64>(j - 1ULL);
                }
            }
        }

        const u64 ans = static_cast<u64>(sum);
        cache.emplace(x, ans);
        return ans;
    }

private:
    using i64 = std::int64_t;

    static constexpr int kSmallPhiPrimes = 7;   // 2,3,5,7,11,13,17
    static constexpr int kSmallPhiMod = 510510;
    static constexpr int kPhiCacheSMax = 128;
    static constexpr u64 kPhiCacheXMax = 50'000'000'000ULL;

    const PrimeTable& table_;
    std::vector<std::array<u32, kSmallPhiPrimes + 1>> phi_table_;

    void build_phi_table() {
        phi_table_.assign(static_cast<std::size_t>(kSmallPhiMod) + 1ULL, {});
        for (int n = 0; n <= kSmallPhiMod; ++n) {
            phi_table_[static_cast<std::size_t>(n)][0] = static_cast<u32>(n);
        }

        for (int s = 1; s <= kSmallPhiPrimes; ++s) {
            const int p = table_.primes[static_cast<std::size_t>(s - 1)];
            for (int n = 0; n <= kSmallPhiMod; ++n) {
                phi_table_[static_cast<std::size_t>(n)][static_cast<std::size_t>(s)] =
                    phi_table_[static_cast<std::size_t>(n)][static_cast<std::size_t>(s - 1)] -
                    phi_table_[static_cast<std::size_t>(n / p)][static_cast<std::size_t>(s - 1)];
            }
        }
    }

    static std::unordered_map<u64, u64>& pi_cache() {
        static thread_local std::unordered_map<u64, u64> cache;
        return cache;
    }

    static std::unordered_map<u64, u64>& phi_cache() {
        static thread_local std::unordered_map<u64, u64> cache;
        return cache;
    }

    u64 phi(u64 x, int s) const {
        if (s == 0) {
            return x;
        }

        if (s <= kSmallPhiPrimes) {
            const u64 block = x / static_cast<u64>(kSmallPhiMod);
            const u64 rem = x % static_cast<u64>(kSmallPhiMod);
            return block * static_cast<u64>(phi_table_[kSmallPhiMod][static_cast<std::size_t>(s)]) +
                   static_cast<u64>(
                       phi_table_[static_cast<std::size_t>(rem)][static_cast<std::size_t>(s)]);
        }

        if (x <= static_cast<u64>(table_.limit) &&
            static_cast<u64>(table_.primes[static_cast<std::size_t>(s - 1)]) *
                    static_cast<u64>(table_.primes[static_cast<std::size_t>(s - 1)]) >
                x) {
            return static_cast<u64>(table_.pi_small[static_cast<std::size_t>(x)] - s + 1);
        }

        const u64 ps = static_cast<u64>(table_.primes[static_cast<std::size_t>(s - 1)]);
        if (ps * ps > x) {
            return pi(x) - static_cast<u64>(s) + 1ULL;
        }

        if (s <= kPhiCacheSMax && x <= kPhiCacheXMax) {
            const u64 key = (x << 8U) ^ static_cast<u64>(s);
            auto& cache = phi_cache();
            const auto it = cache.find(key);
            if (it != cache.end()) {
                return it->second;
            }
            const u64 ans = phi(x, s - 1) - phi(x / ps, s - 1);
            cache.emplace(key, ans);
            return ans;
        }

        return phi(x, s - 1) - phi(x / ps, s - 1);
    }
};

struct FsfMemoKey {
    u16 mask_pos = 0U;
    u64 packed_state = 0ULL;

    bool operator==(const FsfMemoKey& other) const {
        return mask_pos == other.mask_pos && packed_state == other.packed_state;
    }
};

struct FsfMemoKeyHash {
    std::size_t operator()(const FsfMemoKey& key) const {
        const u64 mixed = key.packed_state ^ (static_cast<u64>(key.mask_pos) * 0x9E3779B97F4A7C15ULL);
        return static_cast<std::size_t>(mixed ^ (mixed >> 33U));
    }
};

u64 pack_small_vector(const std::array<u8, 10>& values, int size) {
    u64 packed = 0ULL;
    for (int i = 0; i < size; ++i) {
        packed |= static_cast<u64>(values[static_cast<std::size_t>(i)]) << (6U * i);
    }
    return packed;
}

u64 fsf_from_signature(const std::vector<u8>& signature) {
    const int r = static_cast<int>(signature.size());
    if (r == 0) {
        return 0ULL;
    }

    std::vector<int> masks;
    masks.reserve((1 << r) - 1);
    for (int mask = 1; mask < (1 << r); ++mask) {
        masks.push_back(mask);
    }
    std::sort(masks.begin(), masks.end(), [](int a, int b) {
        const int pa = __builtin_popcount(static_cast<unsigned>(a));
        const int pb = __builtin_popcount(static_cast<unsigned>(b));
        if (pa != pb) {
            return pa > pb;
        }
        return a < b;
    });

    std::vector<std::vector<u8>> mask_bits;
    mask_bits.reserve(masks.size());
    for (const int mask : masks) {
        std::vector<u8> bits;
        for (int i = 0; i < r; ++i) {
            if ((mask >> i) & 1) {
                bits.push_back(static_cast<u8>(i));
            }
        }
        mask_bits.push_back(std::move(bits));
    }

    std::array<u8, 10> remaining{};
    for (int i = 0; i < r; ++i) {
        remaining[static_cast<std::size_t>(i)] = signature[static_cast<std::size_t>(i)];
    }

    std::unordered_map<FsfMemoKey, u64, FsfMemoKeyHash> memo;
    memo.reserve(4096);

    std::function<u64(u16)> dfs = [&](u16 pos) -> u64 {
        if (pos == static_cast<u16>(mask_bits.size())) {
            for (int i = 0; i < r; ++i) {
                if (remaining[static_cast<std::size_t>(i)] != 0U) {
                    return 0ULL;
                }
            }
            return 1ULL;
        }

        const FsfMemoKey key{pos, pack_small_vector(remaining, r)};
        const auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        const std::vector<u8>& bits = mask_bits[static_cast<std::size_t>(pos)];
        u8 max_take = std::numeric_limits<u8>::max();
        for (const u8 bit : bits) {
            max_take = std::min(max_take, remaining[static_cast<std::size_t>(bit)]);
        }

        u64 total = 0ULL;
        for (u32 take = 0U;; ++take) {
            total += dfs(static_cast<u16>(pos + 1U));
            if (take == max_take) {
                break;
            }
            for (const u8 bit : bits) {
                --remaining[static_cast<std::size_t>(bit)];
            }
        }

        for (const u8 bit : bits) {
            remaining[static_cast<std::size_t>(bit)] =
                static_cast<u8>(remaining[static_cast<std::size_t>(bit)] + max_take);
        }

        memo.emplace(key, total);
        return total;
    };

    return dfs(0U);
}

void generate_signatures_dfs(const std::vector<int>& seed_primes,
                             int prime_index,
                             int last_exp,
                             u64 current,
                             u64 limit,
                             std::vector<u8>& cur,
                             std::vector<std::vector<u8>>& out) {
    if (prime_index >= static_cast<int>(seed_primes.size())) {
        return;
    }

    const u64 p = static_cast<u64>(seed_primes[static_cast<std::size_t>(prime_index)]);
    u64 p_pow = 1ULL;

    for (int e = 1; e <= last_exp; ++e) {
        if (p_pow > limit / p) {
            break;
        }
        p_pow *= p;
        if (current > limit / p_pow) {
            break;
        }

        const u64 next_value = current * p_pow;
        cur.push_back(static_cast<u8>(e));
        out.push_back(cur);
        generate_signatures_dfs(seed_primes, prime_index + 1, e, next_value, limit, cur, out);
        cur.pop_back();
    }
}

std::vector<std::vector<u8>> generate_signatures(u64 limit, const PrimeTable& table) {
    std::vector<std::vector<u8>> signatures;
    std::vector<u8> cur;
    generate_signatures_dfs(table.primes, 0, 63, 1ULL, limit, cur, signatures);
    return signatures;
}

void generate_unique_permutations_rec(std::array<int, 64>& freq,
                                      std::vector<u8>& cur,
                                      std::size_t pos,
                                      std::vector<std::vector<u8>>& out) {
    if (pos == cur.size()) {
        out.push_back(cur);
        return;
    }

    for (int value = 1; value < static_cast<int>(freq.size()); ++value) {
        if (freq[static_cast<std::size_t>(value)] == 0) {
            continue;
        }
        --freq[static_cast<std::size_t>(value)];
        cur[pos] = static_cast<u8>(value);
        generate_unique_permutations_rec(freq, cur, pos + 1ULL, out);
        ++freq[static_cast<std::size_t>(value)];
    }
}

std::vector<std::vector<u8>> generate_unique_permutations(const std::vector<u8>& signature) {
    std::array<int, 64> freq{};
    for (const u8 e : signature) {
        ++freq[static_cast<std::size_t>(e)];
    }

    std::vector<std::vector<u8>> out;
    std::vector<u8> cur(signature.size(), 0U);
    generate_unique_permutations_rec(freq, cur, 0ULL, out);
    return out;
}

struct CountMemoKey {
    u64 limit = 0ULL;
    u32 start_index = 0U;
    u8 pos = 0U;

    bool operator==(const CountMemoKey& other) const {
        return limit == other.limit && start_index == other.start_index && pos == other.pos;
    }
};

struct CountMemoKeyHash {
    std::size_t operator()(const CountMemoKey& key) const {
        u64 h = key.limit;
        h ^= static_cast<u64>(key.start_index) * 0x9E3779B97F4A7C15ULL;
        h ^= static_cast<u64>(key.pos) * 0xC2B2AE3D27D4EB4FULL;
        h ^= (h >> 31U);
        return static_cast<std::size_t>(h);
    }
};

class OrderedPrimeTupleCounter {
public:
    OrderedPrimeTupleCounter(const PrimeTable& table, const PrimeCounter& prime_counter)
        : table_(table), prime_counter_(prime_counter) {}

    u64 count_for_exponents(const std::vector<u8>& exponents, u64 limit) const {
        const int r = static_cast<int>(exponents.size());
        if (r == 0) {
            return 0ULL;
        }
        if (r == 1) {
            const u64 root = kth_root_floor(limit, static_cast<int>(exponents[0]));
            return prime_counter_.pi(root);
        }

        std::array<int, 10> tail_sum{};
        int running = 0;
        for (int i = r - 1; i >= 0; --i) {
            running += static_cast<int>(exponents[static_cast<std::size_t>(i)]);
            tail_sum[static_cast<std::size_t>(i)] = running;
        }

        std::unordered_map<CountMemoKey, u64, CountMemoKeyHash> memo;
        memo.reserve(4096);

        std::function<u64(int, u32, u64)> dfs = [&](int pos, u32 start_index, u64 rem_limit) -> u64 {
            if (pos == r) {
                return 1ULL;
            }

            const CountMemoKey key{rem_limit, start_index, static_cast<u8>(pos)};
            const auto it = memo.find(key);
            if (it != memo.end()) {
                return it->second;
            }

            u64 result = 0ULL;

            if (pos == r - 1) {
                const int exp = static_cast<int>(exponents[static_cast<std::size_t>(pos)]);
                const u64 root = kth_root_floor(rem_limit, exp);
                const u64 upto = prime_counter_.pi(root);
                result = (upto > static_cast<u64>(start_index))
                             ? (upto - static_cast<u64>(start_index))
                             : 0ULL;
                memo.emplace(key, result);
                return result;
            }

            const int tail = tail_sum[static_cast<std::size_t>(pos)];
            const u64 max_p = kth_root_floor(rem_limit, tail);
            u64 upto = prime_counter_.pi(max_p);
            if (upto > table_.primes.size()) {
                upto = static_cast<u64>(table_.primes.size());
            }

            const int exp = static_cast<int>(exponents[static_cast<std::size_t>(pos)]);
            for (u64 idx = static_cast<u64>(start_index); idx < upto; ++idx) {
                const u64 p = static_cast<u64>(table_.primes[static_cast<std::size_t>(idx)]);
                const u64 p_pow = pow_u64_capped(p, exp, rem_limit);
                if (p_pow > rem_limit) {
                    break;
                }
                result += dfs(pos + 1, static_cast<u32>(idx + 1ULL), rem_limit / p_pow);
            }

            memo.emplace(key, result);
            return result;
        };

        return dfs(0, 0U, limit);
    }

private:
    const PrimeTable& table_;
    const PrimeCounter& prime_counter_;
};

struct Task {
    std::vector<u8> exponents;
    u64 fsf_value = 0ULL;
};

u128 solve_fast(u64 limit, bool allow_multithreading, unsigned requested_threads) {
    const u64 root_limit = std::max<u64>(isqrt_u64(limit), 100ULL);
    if (root_limit > static_cast<u64>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Sieve limit too large for this implementation.");
    }

    const PrimeTable table(static_cast<int>(root_limit));
    const PrimeCounter prime_counter(table);
    const OrderedPrimeTupleCounter tuple_counter(table, prime_counter);

    const std::vector<std::vector<u8>> signatures = generate_signatures(limit, table);

    std::vector<Task> tasks;
    tasks.reserve(16'000ULL);
    for (const auto& signature : signatures) {
        const u64 fsf = fsf_from_signature(signature);
        const std::vector<std::vector<u8>> permutations = generate_unique_permutations(signature);
        for (const auto& p : permutations) {
            tasks.push_back(Task{p, fsf});
        }
    }

    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, tasks.size());

    if (threads == 1U) {
        u128 total = 0;
        for (const Task& task : tasks) {
            const u64 count = tuple_counter.count_for_exponents(task.exponents, limit);
            total += static_cast<u128>(task.fsf_value) * static_cast<u128>(count);
        }
        return total;
    }

    std::atomic<std::size_t> next_task{0ULL};
    std::vector<u128> partials(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            const OrderedPrimeTupleCounter local_counter(table, prime_counter);
            u128 local = 0;
            while (true) {
                const std::size_t idx = next_task.fetch_add(1ULL, std::memory_order_relaxed);
                if (idx >= tasks.size()) {
                    break;
                }
                const Task& task = tasks[idx];
                const u64 count = local_counter.count_for_exponents(task.exponents, limit);
                local += static_cast<u128>(task.fsf_value) * static_cast<u128>(count);
            }
            partials[static_cast<std::size_t>(tid)] = local;
        });
    }

    for (auto& th : workers) {
        th.join();
    }

    u128 total = 0;
    for (const u128 value : partials) {
        total += value;
    }
    return total;
}

std::vector<int> squarefree_numbers(int limit) {
    std::vector<int> sqf;
    std::vector<int> square_div(static_cast<std::size_t>(limit) + 1ULL, 0);
    for (int p = 2; p * p <= limit; ++p) {
        const int pp = p * p;
        for (int n = pp; n <= limit; n += pp) {
            square_div[static_cast<std::size_t>(n)] = 1;
        }
    }
    for (int n = 2; n <= limit; ++n) {
        if (square_div[static_cast<std::size_t>(n)] == 0) {
            sqf.push_back(n);
        }
    }
    return sqf;
}

u64 brute_fsf_single(u64 n,
                     std::size_t start_idx,
                     const std::vector<int>& sqf,
                     std::unordered_map<u64, u64>& memo) {
    const u64 key = (n << 16U) ^ static_cast<u64>(start_idx);
    const auto it = memo.find(key);
    if (it != memo.end()) {
        return it->second;
    }

    u64 total = 0ULL;
    for (std::size_t i = start_idx; i < sqf.size(); ++i) {
        const u64 d = static_cast<u64>(sqf[i]);
        if (d > n) {
            break;
        }
        if (n % d != 0ULL) {
            continue;
        }
        if (d == n) {
            ++total;
        } else {
            total += brute_fsf_single(n / d, i, sqf, memo);
        }
    }

    memo.emplace(key, total);
    return total;
}

u64 brute_s(u64 limit) {
    if (limit < 2ULL) {
        return 0ULL;
    }
    if (limit > 100'000ULL) {
        throw std::runtime_error("brute_s is intended only for small limits.");
    }

    const std::vector<int> sqf = squarefree_numbers(static_cast<int>(limit));
    u64 total = 0ULL;
    for (u64 n = 2ULL; n <= limit; ++n) {
        std::unordered_map<u64, u64> memo;
        total += brute_fsf_single(n, 0ULL, sqf, memo);
    }
    return total;
}

void run_checkpoints() {
    const u64 brute_100 = brute_s(kCheckpointKnownLimit);
    if (brute_100 != kCheckpointKnownExpected) {
        throw std::runtime_error("Brute-force checkpoint S(100)=193 failed.");
    }

    const u64 fast_100 = static_cast<u64>(solve_fast(kCheckpointKnownLimit, false, 1U));
    if (fast_100 != kCheckpointKnownExpected) {
        throw std::runtime_error("Fast checkpoint S(100)=193 failed.");
    }

    const u64 brute_small = brute_s(kCheckpointBruteLimit);
    const u64 fast_small = static_cast<u64>(solve_fast(kCheckpointBruteLimit, false, 1U));
    if (fast_small != brute_small) {
        throw std::runtime_error("Brute-force consistency checkpoint failed.");
    }

    const u64 single = static_cast<u64>(solve_fast(kThreadConsistencyLimit, false, 1U));
    const u64 threaded = static_cast<u64>(solve_fast(kThreadConsistencyLimit, true, 2U));
    if (single != threaded) {
        throw std::runtime_error("Thread consistency checkpoint failed.");
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    try {
        if (options.run_checkpoints) {
            run_checkpoints();
        }

        const auto start = std::chrono::steady_clock::now();
        const u128 answer = solve_fast(
            options.limit, options.allow_multithreading, options.requested_threads);
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - start;

        std::cout << "S(" << options.limit << ") = " << to_string_u128(answer) << '\n';
        std::cout << "Elapsed: " << elapsed.count() << " s\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
