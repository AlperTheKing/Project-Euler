#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

struct Options {
    int limit = 80;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

struct PrimeTopGroup {
    int prime = 0;
    int prime_power = 1;
    std::vector<int> numbers;
    std::vector<u64> valid_masks;
    u64 used_union_mask = 0ULL;
};

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

    std::uint64_t parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<unsigned>::max())) {
            return false;
        }
    }

    value = static_cast<unsigned>(parsed);
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

    long long parsed = 0LL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10LL + static_cast<long long>(c - '0');
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
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
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

std::vector<int> sieve_primes(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit) + 1ULL, true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }

    for (int p = 2; p * p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            primes.push_back(p);
        }
    }
    return primes;
}

i64 checked_lcm(const i64 a, const i64 b) {
    const i64 g = std::gcd(a, b);
    const __int128 scaled = static_cast<__int128>(a / g) * static_cast<__int128>(b);
    if (scaled > static_cast<__int128>(std::numeric_limits<i64>::max())) {
        throw std::overflow_error("LCM overflow");
    }
    return static_cast<i64>(scaled);
}

i64 mod_inverse(i64 a, i64 mod) {
    i64 t = 0;
    i64 new_t = 1;
    i64 r = mod;
    i64 new_r = a % mod;

    while (new_r != 0) {
        const i64 q = r / new_r;

        const i64 next_t = t - q * new_t;
        t = new_t;
        new_t = next_t;

        const i64 next_r = r - q * new_r;
        r = new_r;
        new_r = next_r;
    }

    if (r != 1) {
        throw std::runtime_error("Modular inverse does not exist");
    }

    t %= mod;
    if (t < 0) {
        t += mod;
    }
    return t;
}

std::vector<PrimeTopGroup> build_odd_prime_top_groups(const int limit) {
    const std::vector<int> primes = sieve_primes(limit);
    std::vector<PrimeTopGroup> groups;

    for (const int p : primes) {
        if (p == 2) {
            continue;
        }

        int prime_power = p;
        while (prime_power <= limit / p) {
            prime_power *= p;
        }

        PrimeTopGroup group;
        group.prime = p;
        group.prime_power = prime_power;

        const int next_power = prime_power * p;
        for (int n = 2; n <= limit; ++n) {
            if (n % prime_power != 0) {
                continue;
            }
            if (next_power <= limit && n % next_power == 0) {
                continue;
            }
            group.numbers.push_back(n);
        }

        if (group.numbers.empty()) {
            continue;
        }

        if (group.numbers.size() > 20U) {
            throw std::runtime_error("Unexpectedly large odd-prime top group");
        }

        // For odd prime p, modulo p^2 only numbers with maximal p-adic valuation survive.
        // Their normalized coefficients are inverses of (n / p^e)^2 modulo p^2.
        const i64 mod = static_cast<i64>(p) * static_cast<i64>(p);
        std::vector<i64> coeffs;
        coeffs.reserve(group.numbers.size());
        for (const int n : group.numbers) {
            const i64 m = static_cast<i64>(n / prime_power);
            const i64 sq = (m * m) % mod;
            coeffs.push_back(mod_inverse(sq, mod));
        }

        const u64 mask_count = 1ULL << group.numbers.size();
        group.valid_masks.reserve(static_cast<std::size_t>(mask_count));

        for (u64 mask = 0ULL; mask < mask_count; ++mask) {
            i64 residue = 0;
            for (std::size_t i = 0; i < coeffs.size(); ++i) {
                if ((mask >> i) & 1ULL) {
                    residue += coeffs[i];
                    residue %= mod;
                }
            }
            if (residue == 0) {
                group.valid_masks.push_back(mask);
                group.used_union_mask |= mask;
            }
        }

        groups.push_back(std::move(group));
    }

    return groups;
}

std::unordered_map<i64, u64> subset_sum_counts(const std::vector<i64>& weights) {
    if (weights.size() > 30U) {
        throw std::runtime_error("Subset half is unexpectedly large");
    }

    const u64 mask_count = 1ULL << weights.size();
    std::unordered_map<i64, u64> counts;
    counts.reserve(static_cast<std::size_t>(mask_count * 1.3));

    for (u64 mask = 0ULL; mask < mask_count; ++mask) {
        i64 sum = 0;
        for (std::size_t i = 0; i < weights.size(); ++i) {
            if ((mask >> i) & 1ULL) {
                sum += weights[i];
            }
        }
        ++counts[sum];
    }

    return counts;
}

u64 count_with_meet_in_middle(const std::vector<i64>& free_weights,
                              const std::vector<std::pair<i64, u64>>& base_states,
                              const i64 target,
                              const unsigned thread_count) {
    const std::size_t mid = free_weights.size() / 2U;
    const std::vector<i64> left_weights(free_weights.begin(), free_weights.begin() + static_cast<std::ptrdiff_t>(mid));
    const std::vector<i64> right_weights(free_weights.begin() + static_cast<std::ptrdiff_t>(mid), free_weights.end());

    const std::unordered_map<i64, u64> left_counts = subset_sum_counts(left_weights);
    const std::unordered_map<i64, u64> right_counts = subset_sum_counts(right_weights);

    std::vector<std::pair<i64, u64>> left_entries;
    left_entries.reserve(left_counts.size());
    for (const auto& [sum, ways] : left_counts) {
        left_entries.push_back({sum, ways});
    }

    auto count_chunk = [&](const std::size_t begin, const std::size_t end) {
        u64 local = 0ULL;
        for (std::size_t i = begin; i < end; ++i) {
            const i64 left_sum = left_entries[i].first;
            const u64 left_ways = left_entries[i].second;

            for (const auto& [base_sum, base_ways] : base_states) {
                const i64 need = target - base_sum - left_sum;
                const auto it = right_counts.find(need);
                if (it == right_counts.end()) {
                    continue;
                }
                local += left_ways * it->second * base_ways;
            }
        }
        return local;
    };

    if (thread_count <= 1U || left_entries.size() < 12000U) {
        return count_chunk(0U, left_entries.size());
    }

    const unsigned workers = std::min<unsigned>(thread_count, static_cast<unsigned>(left_entries.size()));
    std::atomic<std::size_t> next(0U);
    constexpr std::size_t kChunk = 256U;

    std::vector<u64> partial(workers, 0ULL);
    std::vector<std::thread> pool;
    pool.reserve(workers);

    for (unsigned t = 0U; t < workers; ++t) {
        pool.emplace_back([&, t]() {
            u64 local = 0ULL;
            while (true) {
                const std::size_t begin = next.fetch_add(kChunk, std::memory_order_relaxed);
                if (begin >= left_entries.size()) {
                    break;
                }
                const std::size_t end = std::min(begin + kChunk, left_entries.size());
                local += count_chunk(begin, end);
            }
            partial[t] = local;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u64 total = 0ULL;
    for (const u64 v : partial) {
        total += v;
    }
    return total;
}

u64 count_representations(const int limit, const unsigned thread_count) {
    if (limit < 2) {
        return 0ULL;
    }
    if (limit > 80) {
        throw std::runtime_error("This implementation targets limits up to 80.");
    }

    const std::vector<PrimeTopGroup> groups = build_odd_prime_top_groups(limit);

    std::vector<int> number_to_group(static_cast<std::size_t>(limit) + 1ULL, -1);
    std::vector<int> number_to_group_bit(static_cast<std::size_t>(limit) + 1ULL, -1);

    for (std::size_t g = 0; g < groups.size(); ++g) {
        for (std::size_t i = 0; i < groups[g].numbers.size(); ++i) {
            const int n = groups[g].numbers[i];
            if (number_to_group[static_cast<std::size_t>(n)] != -1) {
                throw std::runtime_error("Odd-prime top groups overlapped unexpectedly.");
            }
            number_to_group[static_cast<std::size_t>(n)] = static_cast<int>(g);
            number_to_group_bit[static_cast<std::size_t>(n)] = static_cast<int>(i);
        }
    }

    std::vector<const PrimeTopGroup*> active_groups;
    active_groups.reserve(groups.size());

    for (const PrimeTopGroup& group : groups) {
        if (group.valid_masks.size() == 1U && group.valid_masks.front() == 0ULL) {
            continue;
        }
        active_groups.push_back(&group);
    }

    std::vector<char> number_is_possible(static_cast<std::size_t>(limit) + 1ULL, false);
    std::vector<int> free_numbers;
    free_numbers.reserve(static_cast<std::size_t>(limit));

    for (int n = 2; n <= limit; ++n) {
        const int gid = number_to_group[static_cast<std::size_t>(n)];
        if (gid == -1) {
            number_is_possible[static_cast<std::size_t>(n)] = true;
            free_numbers.push_back(n);
            continue;
        }

        const PrimeTopGroup& group = groups[static_cast<std::size_t>(gid)];
        const int bit = number_to_group_bit[static_cast<std::size_t>(n)];
        if (((group.used_union_mask >> bit) & 1ULL) != 0ULL) {
            number_is_possible[static_cast<std::size_t>(n)] = true;
        }
    }

    i64 common = 1;
    for (int n = 2; n <= limit; ++n) {
        if (!number_is_possible[static_cast<std::size_t>(n)]) {
            continue;
        }
        common = checked_lcm(common, static_cast<i64>(n));
    }

    const __int128 sq = static_cast<__int128>(common) * static_cast<__int128>(common);
    if (sq > static_cast<__int128>(std::numeric_limits<i64>::max())) {
        throw std::overflow_error("Common denominator square overflow");
    }
    const i64 common_sq = static_cast<i64>(sq);
    if ((common_sq % 2LL) != 0LL) {
        throw std::runtime_error("Unexpected odd denominator square");
    }
    const i64 target = common_sq / 2LL;

    std::vector<i64> weight_by_number(static_cast<std::size_t>(limit) + 1ULL, 0LL);
    for (int n = 2; n <= limit; ++n) {
        if (!number_is_possible[static_cast<std::size_t>(n)]) {
            continue;
        }
        const i64 nn = static_cast<i64>(n) * static_cast<i64>(n);
        if ((common_sq % nn) != 0LL) {
            throw std::runtime_error("Weight is not integral");
        }
        weight_by_number[static_cast<std::size_t>(n)] = common_sq / nn;
    }

    std::unordered_map<i64, u64> base_states;
    base_states.reserve(8U);
    base_states[0LL] = 1ULL;

    for (const PrimeTopGroup* group_ptr : active_groups) {
        const PrimeTopGroup& group = *group_ptr;

        std::unordered_map<i64, u64> option_counts;
        option_counts.reserve(group.valid_masks.size() * 2ULL + 1ULL);

        for (const u64 mask : group.valid_masks) {
            i64 sum = 0LL;
            for (std::size_t i = 0; i < group.numbers.size(); ++i) {
                if (((mask >> i) & 1ULL) == 0ULL) {
                    continue;
                }
                const int n = group.numbers[i];
                sum += weight_by_number[static_cast<std::size_t>(n)];
            }
            ++option_counts[sum];
        }

        std::unordered_map<i64, u64> next_states;
        next_states.reserve(base_states.size() * option_counts.size() * 2ULL + 1ULL);

        for (const auto& [base_sum, base_ways] : base_states) {
            for (const auto& [opt_sum, opt_ways] : option_counts) {
                next_states[base_sum + opt_sum] += base_ways * opt_ways;
            }
        }

        base_states.swap(next_states);
    }

    std::vector<i64> free_weights;
    free_weights.reserve(free_numbers.size());
    for (const int n : free_numbers) {
        free_weights.push_back(weight_by_number[static_cast<std::size_t>(n)]);
    }

    std::vector<std::pair<i64, u64>> base_entries;
    base_entries.reserve(base_states.size());
    for (const auto& [sum, ways] : base_states) {
        base_entries.push_back({sum, ways});
    }

    return count_with_meet_in_middle(free_weights, base_entries, target, thread_count);
}

u64 brute_force_count(const int limit) {
    if (limit < 2) {
        return 0ULL;
    }

    i64 common = 1;
    for (int n = 2; n <= limit; ++n) {
        common = checked_lcm(common, static_cast<i64>(n));
    }

    const __int128 sq = static_cast<__int128>(common) * static_cast<__int128>(common);
    if (sq > static_cast<__int128>(std::numeric_limits<i64>::max())) {
        throw std::overflow_error("Brute-force denominator square overflow");
    }
    const i64 common_sq = static_cast<i64>(sq);
    if ((common_sq % 2LL) != 0LL) {
        throw std::runtime_error("Unexpected odd brute-force denominator square");
    }
    const i64 target = common_sq / 2LL;

    std::vector<i64> weights;
    weights.reserve(static_cast<std::size_t>(limit - 1));
    for (int n = 2; n <= limit; ++n) {
        weights.push_back(common_sq / (static_cast<i64>(n) * static_cast<i64>(n)));
    }

    if (weights.size() > 26U) {
        throw std::runtime_error("Brute-force checkpoint requested above safe size.");
    }

    const u64 mask_count = 1ULL << weights.size();
    u64 count = 0ULL;
    for (u64 mask = 0ULL; mask < mask_count; ++mask) {
        i64 sum = 0LL;
        for (std::size_t i = 0; i < weights.size(); ++i) {
            if ((mask >> i) & 1ULL) {
                sum += weights[i];
            }
        }
        if (sum == target) {
            ++count;
        }
    }

    return count;
}

bool subset_sums_to_half(const std::vector<int>& subset) {
    i64 common = 2;
    for (const int n : subset) {
        common = checked_lcm(common, static_cast<i64>(n));
    }

    const __int128 sq = static_cast<__int128>(common) * static_cast<__int128>(common);
    if (sq > static_cast<__int128>(std::numeric_limits<i64>::max())) {
        throw std::overflow_error("Checkpoint denominator square overflow");
    }
    const i64 common_sq = static_cast<i64>(sq);

    i64 lhs = 0LL;
    for (const int n : subset) {
        lhs += common_sq / (static_cast<i64>(n) * static_cast<i64>(n));
    }

    return lhs * 2LL == common_sq;
}

bool run_checkpoints(const unsigned thread_count) {
    const std::vector<int> known_one = {2, 3, 4, 5, 7, 12, 15, 20, 28, 35};
    const std::vector<int> known_two = {2, 3, 4, 6, 7, 9, 10, 20, 28, 35, 36, 45};
    const std::vector<int> known_three = {2, 3, 4, 6, 7, 9, 12, 15, 28, 30, 35, 36, 45};

    if (!subset_sums_to_half(known_one)) {
        std::cerr << "Checkpoint failed: first known decomposition is invalid.\n";
        return false;
    }
    if (!subset_sums_to_half(known_two)) {
        std::cerr << "Checkpoint failed: second known decomposition is invalid.\n";
        return false;
    }
    if (!subset_sums_to_half(known_three)) {
        std::cerr << "Checkpoint failed: third known decomposition is invalid.\n";
        return false;
    }

    const u64 count_45 = count_representations(45, thread_count);
    if (count_45 != 3ULL) {
        std::cerr << "Checkpoint failed: expected 3 solutions for 2 <= n <= 45, got "
                  << count_45 << ".\n";
        return false;
    }

    const int brute_limit = 18;
    const u64 brute = brute_force_count(brute_limit);
    const u64 fast = count_representations(brute_limit, 1U);
    if (brute != fast) {
        std::cerr << "Checkpoint failed: fast/brute mismatch at limit " << brute_limit
                  << " (fast=" << fast << ", brute=" << brute << ").\n";
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

    if (options.limit < 2 || options.limit > 80) {
        std::cerr << "Please use --limit in [2, 80].\n";
        return 1;
    }

    const unsigned thread_count = pick_thread_count(options.requested_threads);

    try {
        if (options.run_checkpoints && !run_checkpoints(thread_count)) {
            return 1;
        }

        const u64 answer = count_representations(options.limit, thread_count);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
