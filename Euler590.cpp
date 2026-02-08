#include <algorithm>
#include <array>
#include <cstdint>
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

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 MOD = 1'000'000'000U;
constexpr u32 POW2_THRESHOLD = 9U;
constexpr u32 POW2_CYCLE = 1'562'500U;  // ord_{5^9}(2) = 4 * 5^8

struct Options {
    int n = 50'000;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

struct ChoiceData {
    u32 coef = 0U;      // (-1)^i * C(count, i) mod MOD
    u32 exp_mod = 0U;   // factor exponent mod POW2_CYCLE
    u8 exp_cap = 0U;    // min(real exponent factor, 9)
};

struct GroupData {
    int d = 0;          // d = e_p + 1
    int count = 0;      // how many primes have this d
    std::vector<ChoiceData> choices;
};

struct RestState {
    u32 coef = 0U;
    u32 exp_mod = 0U;
    u8 exp_cap = 1U;
};

struct Pow2Tables {
    std::array<u32, POW2_THRESHOLD> small{};
    std::vector<u32> cycle_shifted;  // cycle_shifted[r] = 2^(9+r) mod MOD
};

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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
        if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
            return false;
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
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

u8 cap_mul_9(const u8 a, const u8 b) {
    const unsigned product = static_cast<unsigned>(a) * static_cast<unsigned>(b);
    if (product >= POW2_THRESHOLD) {
        return static_cast<u8>(POW2_THRESHOLD);
    }
    return static_cast<u8>(product);
}

u32 pow_mod_u32(u32 base, u64 exp, const u32 mod) {
    if (mod == 1U) {
        return 0U;
    }

    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = (result * cur) % mod;
        }
        cur = (cur * cur) % mod;
        exp >>= 1ULL;
    }
    return static_cast<u32>(result);
}

const Pow2Tables& get_pow2_tables() {
    static const Pow2Tables tables = []() {
        Pow2Tables t;
        t.small[0] = 1U;
        for (u32 i = 1U; i < POW2_THRESHOLD; ++i) {
            t.small[i] = (2U * t.small[i - 1U]) % MOD;
        }

        t.cycle_shifted.assign(POW2_CYCLE, 0U);
        u32 value = static_cast<u32>((2ULL * t.small[POW2_THRESHOLD - 1U]) % MOD);  // 2^9 mod MOD
        for (u32 i = 0U; i < POW2_CYCLE; ++i) {
            t.cycle_shifted[i] = value;
            value = static_cast<u32>((2ULL * value) % MOD);
        }
        return t;
    }();
    return tables;
}

u32 pow2_mod_from_exp(const u32 exp_mod, const u8 exp_cap) {
    const Pow2Tables& tables = get_pow2_tables();
    if (exp_cap < POW2_THRESHOLD) {
        return tables.small[exp_cap];
    }

    // For E >= 9:
    // 2^E mod 10^9 = 2^(9 + ((E - 9) mod ord_{5^9}(2))) mod 10^9.
    u32 idx = exp_mod;
    if (idx >= POW2_THRESHOLD) {
        idx -= POW2_THRESHOLD;
    } else {
        idx += POW2_CYCLE - POW2_THRESHOLD;
    }
    return tables.cycle_shifted[idx];
}

std::vector<int> sieve_primes(const int limit) {
    if (limit < 2) {
        return {};
    }

    std::vector<bool> is_prime(static_cast<std::size_t>(limit) + 1ULL, true);
    is_prime[0] = false;
    is_prime[1] = false;

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

std::vector<int> exponent_plus_one_list_for_lcm_upto_n(const int n) {
    const std::vector<int> primes = sieve_primes(n);
    std::vector<int> result;
    result.reserve(primes.size());

    for (const int p : primes) {
        int d = 0;
        int x = n;
        while (x > 0) {
            x /= p;
            ++d;
        }
        result.push_back(d);
    }

    return result;
}

std::vector<u32> build_binomial_row_mod(const int n) {
    std::vector<u32> row(static_cast<std::size_t>(n) + 1ULL, 0U);
    row[0] = 1U;

    for (int i = 1; i <= n; ++i) {
        for (int k = i; k >= 1; --k) {
            u32 value = row[static_cast<std::size_t>(k)] +
                        row[static_cast<std::size_t>(k - 1)];
            if (value >= MOD) {
                value -= MOD;
            }
            row[static_cast<std::size_t>(k)] = value;
        }
    }

    return row;
}

GroupData build_group_data(const int d, const int count) {
    GroupData group;
    group.d = d;
    group.count = count;
    group.choices.assign(static_cast<std::size_t>(count) + 1ULL, ChoiceData{});

    const std::vector<u32> binom = build_binomial_row_mod(count);

    std::vector<u32> pow_d_mod(static_cast<std::size_t>(count) + 1ULL, 1U);
    std::vector<u32> pow_dm1_mod(static_cast<std::size_t>(count) + 1ULL, 1U);
    std::vector<u8> pow_d_cap(static_cast<std::size_t>(count) + 1ULL, 1U);
    std::vector<u8> pow_dm1_cap(static_cast<std::size_t>(count) + 1ULL, 1U);

    for (int k = 1; k <= count; ++k) {
        pow_d_mod[static_cast<std::size_t>(k)] = static_cast<u32>(
            (static_cast<u64>(pow_d_mod[static_cast<std::size_t>(k - 1)]) *
             static_cast<u64>(d)) %
            POW2_CYCLE);
        pow_dm1_mod[static_cast<std::size_t>(k)] = static_cast<u32>(
            (static_cast<u64>(pow_dm1_mod[static_cast<std::size_t>(k - 1)]) *
             static_cast<u64>(d - 1)) %
            POW2_CYCLE);

        pow_d_cap[static_cast<std::size_t>(k)] = cap_mul_9(
            pow_d_cap[static_cast<std::size_t>(k - 1)], static_cast<u8>(d));
        pow_dm1_cap[static_cast<std::size_t>(k)] = cap_mul_9(
            pow_dm1_cap[static_cast<std::size_t>(k - 1)], static_cast<u8>(d - 1));
    }

    for (int i = 0; i <= count; ++i) {
        ChoiceData choice;

        const u32 c = binom[static_cast<std::size_t>(i)];
        if ((i & 1) == 0) {
            choice.coef = c;
        } else {
            choice.coef = (c == 0U ? 0U : MOD - c);
        }

        const u32 left = pow_dm1_mod[static_cast<std::size_t>(i)];
        const u32 right = pow_d_mod[static_cast<std::size_t>(count - i)];
        choice.exp_mod = static_cast<u32>(
            (static_cast<u64>(left) * static_cast<u64>(right)) % POW2_CYCLE);

        choice.exp_cap = cap_mul_9(pow_dm1_cap[static_cast<std::size_t>(i)],
                                   pow_d_cap[static_cast<std::size_t>(count - i)]);

        group.choices[static_cast<std::size_t>(i)] = choice;
    }

    return group;
}

std::vector<GroupData> build_groups_for_hl(const int n) {
    const std::vector<int> d_list = exponent_plus_one_list_for_lcm_upto_n(n);

    std::unordered_map<int, int> count_by_d;
    count_by_d.reserve(d_list.size() * 2ULL + 1ULL);
    for (const int d : d_list) {
        ++count_by_d[d];
    }

    std::vector<std::pair<int, int>> sorted_pairs;
    sorted_pairs.reserve(count_by_d.size());
    for (const auto& kv : count_by_d) {
        sorted_pairs.push_back(kv);
    }

    // Sort by multiplicity descending so the largest bucket is parallelized.
    std::sort(sorted_pairs.begin(), sorted_pairs.end(),
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                  if (a.second != b.second) {
                      return a.second > b.second;
                  }
                  return a.first > b.first;
              });

    std::vector<GroupData> groups;
    groups.reserve(sorted_pairs.size());
    for (const auto& [d, count] : sorted_pairs) {
        groups.push_back(build_group_data(d, count));
    }

    return groups;
}

std::vector<RestState> build_rest_states(const std::vector<GroupData>& groups) {
    std::vector<RestState> states;
    states.push_back(RestState{1U, 1U, 1U});

    for (std::size_t g = 1; g < groups.size(); ++g) {
        const GroupData& group = groups[g];
        std::vector<RestState> next;
        next.reserve(states.size() * group.choices.size());

        for (const RestState& state : states) {
            for (const ChoiceData& choice : group.choices) {
                const u32 coef = static_cast<u32>(
                    (static_cast<u64>(state.coef) * static_cast<u64>(choice.coef)) %
                    MOD);
                if (coef == 0U) {
                    continue;
                }

                const u32 exp_mod = static_cast<u32>(
                    (static_cast<u64>(state.exp_mod) *
                     static_cast<u64>(choice.exp_mod)) %
                    POW2_CYCLE);
                const u8 exp_cap = cap_mul_9(state.exp_cap, choice.exp_cap);
                next.push_back(RestState{coef, exp_mod, exp_cap});
            }
        }

        states.swap(next);
    }

    return states;
}

u32 solve_hl_mod(const int n, const unsigned requested_threads) {
    const std::vector<GroupData> groups = build_groups_for_hl(n);
    if (groups.empty()) {
        return 0U;
    }

    const GroupData& lead = groups[0];
    const std::vector<RestState> rest_states = build_rest_states(groups);

    const std::size_t lead_choices = lead.choices.size();
    unsigned thread_count = requested_threads;
    if (thread_count == 0U) {
        thread_count = 1U;
    }
    thread_count = std::min<unsigned>(thread_count, static_cast<unsigned>(lead_choices));
    if (thread_count == 0U) {
        thread_count = 1U;
    }

    std::vector<u64> partial(static_cast<std::size_t>(thread_count), 0ULL);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    const std::size_t block = (lead_choices + static_cast<std::size_t>(thread_count) - 1ULL) /
                              static_cast<std::size_t>(thread_count);

    auto worker = [&](const unsigned worker_id) {
        const std::size_t begin = static_cast<std::size_t>(worker_id) * block;
        const std::size_t end = std::min<std::size_t>(lead_choices, begin + block);
        u64 local = 0ULL;

        for (std::size_t i = begin; i < end; ++i) {
            const ChoiceData& prefix = lead.choices[i];
            for (const RestState& state : rest_states) {
                const u32 coef = static_cast<u32>(
                    (static_cast<u64>(prefix.coef) * static_cast<u64>(state.coef)) %
                    MOD);
                if (coef == 0U) {
                    continue;
                }

                const u8 exp_cap = cap_mul_9(prefix.exp_cap, state.exp_cap);
                const u32 exp_mod = static_cast<u32>(
                    (static_cast<u64>(prefix.exp_mod) * static_cast<u64>(state.exp_mod)) %
                    POW2_CYCLE);
                const u32 value_2exp = pow2_mod_from_exp(exp_mod, exp_cap);

                local += (static_cast<u64>(coef) * static_cast<u64>(value_2exp)) % MOD;
            }
        }

        partial[static_cast<std::size_t>(worker_id)] = local % MOD;
    };

    for (unsigned t = 0U; t < thread_count; ++t) {
        workers.emplace_back(worker, t);
    }
    for (std::thread& th : workers) {
        th.join();
    }

    u64 total = 0ULL;
    for (const u64 value : partial) {
        total += value;
    }
    return static_cast<u32>(total % MOD);
}

u64 brute_h(const u64 n) {
    std::vector<u64> divisors;
    for (u64 d = 1ULL; d <= n; ++d) {
        if (n % d == 0ULL) {
            divisors.push_back(d);
        }
    }

    const std::size_t m = divisors.size();
    if (m >= 63U) {
        throw std::runtime_error("brute_h requested for too many divisors");
    }

    const u64 subset_count = (1ULL << m);
    u64 answer = 0ULL;
    for (u64 mask = 1ULL; mask < subset_count; ++mask) {
        u64 current_lcm = 1ULL;
        for (std::size_t i = 0; i < m; ++i) {
            if (((mask >> i) & 1ULL) == 0ULL) {
                continue;
            }
            const u64 x = divisors[i];
            const u64 g = std::gcd(current_lcm, x);
            current_lcm = (current_lcm / g) * x;
        }
        if (current_lcm == n) {
            ++answer;
        }
    }

    return answer;
}

u32 solve_hl_direct_ie_small(const int n) {
    const std::vector<int> d_list = exponent_plus_one_list_for_lcm_upto_n(n);
    const int k = static_cast<int>(d_list.size());
    if (k >= 60) {
        throw std::runtime_error("Direct IE checkpoint is only for small k");
    }

    const u64 subset_count = (1ULL << static_cast<u64>(k));
    u64 answer = 0ULL;

    for (u64 mask = 0ULL; mask < subset_count; ++mask) {
        u64 exponent = 1ULL;
        for (int i = 0; i < k; ++i) {
            const bool in_subset = ((mask >> static_cast<u64>(i)) & 1ULL) != 0ULL;
            const int factor = in_subset ? (d_list[static_cast<std::size_t>(i)] - 1)
                                         : d_list[static_cast<std::size_t>(i)];
            exponent *= static_cast<u64>(factor);
        }

        const u32 term = pow_mod_u32(2U, exponent, MOD);
        if ((__builtin_popcountll(mask) & 1) != 0) {
            answer += (MOD - term);
        } else {
            answer += term;
        }
        answer %= MOD;
    }

    return static_cast<u32>(answer);
}

bool run_checkpoints(const unsigned thread_count) {
    if (brute_h(6ULL) != 10ULL) {
        std::cerr << "Checkpoint failed: H(6) must be 10\n";
        return false;
    }

    const u32 sample = solve_hl_mod(4, 1U);
    if (sample != 44U) {
        std::cerr << "Checkpoint failed: HL(4) must be 44, got " << sample << '\n';
        return false;
    }

    for (int n = 2; n <= 20; ++n) {
        const u32 fast = solve_hl_mod(n, 1U);
        const u32 direct = solve_hl_direct_ie_small(n);
        if (fast != direct) {
            std::cerr << "Checkpoint failed at n=" << n << ": fast=" << fast
                      << " direct=" << direct << '\n';
            return false;
        }
    }

    if (thread_count > 1U) {
        const unsigned check_threads = std::min<unsigned>(thread_count, 4U);
        const u32 single = solve_hl_mod(200, 1U);
        const u32 multi = solve_hl_mod(200, check_threads);
        if (single != multi) {
            std::cerr << "Checkpoint failed: thread consistency mismatch at n=200\n";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.n < 2) {
        std::cerr << "n must be at least 2\n";
        return 1;
    }

    const unsigned thread_count = pick_thread_count(options.requested_threads);

    if (options.run_checkpoints) {
        if (!run_checkpoints(thread_count)) {
            return 1;
        }
    }

    const u32 answer = solve_hl_mod(options.n, thread_count);
    std::cout << answer << '\n';
    return 0;
}
