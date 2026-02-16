#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include <cmath>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = unsigned __int128;

struct Options {
    int target = 150000;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 1;
}

struct U128Hash {
    std::size_t operator()(const u128 v) const noexcept {
        const std::uint64_t lo = static_cast<std::uint64_t>(v);
        const std::uint64_t hi = static_cast<std::uint64_t>(v >> 64);
        return static_cast<std::size_t>(lo ^ (hi * 0x9e3779b97f4a7c15ULL));
    }
};

struct PrimeTable {
    std::vector<u32> primes{2U};
    u32 next_candidate = 3U;

    void ensure_up_to(const u64 limit) {
        while (primes.back() < limit) {
            bool is_prime = true;
            for (const u32 p : primes) {
                const u64 pu = static_cast<u64>(p);
                if (pu * pu > static_cast<u64>(next_candidate)) {
                    break;
                }
                if (next_candidate % p == 0U) {
                    is_prime = false;
                    break;
                }
            }
            if (is_prime) {
                primes.push_back(next_candidate);
            }
            next_candidate += 2U;
        }
    }
};

void build_divisors_leq(
    const std::vector<std::pair<u64, int>>& factors,
    const std::size_t idx,
    const u64 current,
    const u64 limit,
    std::vector<u64>& out
) {
    if (idx == factors.size()) {
        out.push_back(current);
        return;
    }

    const u64 prime = factors[idx].first;
    const int exponent = factors[idx].second;
    u64 power = 1ULL;

    for (int e = 0; e <= exponent; ++e) {
        if (current > limit / power) {
            break;
        }
        build_divisors_leq(factors, idx + 1U, current * power, limit, out);
        if (e == exponent || power > limit / prime) {
            break;
        }
        power *= prime;
    }
}

u64 nth_alexandrian(const int target) {
    std::vector<u128> heap;
    heap.reserve(static_cast<std::size_t>(target) + 16ULL);
    std::unordered_set<u128, U128Hash> active;
    active.reserve(static_cast<std::size_t>(target) * 2ULL + 64ULL);

    PrimeTable table;
    std::vector<std::pair<u64, int>> factors;
    factors.reserve(8U);
    std::vector<u64> divisors;
    divisors.reserve(64U);

    for (u64 p = 1ULL;; ++p) {
        const u64 n = p * p + 1ULL;
        const u64 root = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
        table.ensure_up_to(root + 1ULL);

        factors.clear();
        u64 rem = n;
        for (const u32 prime : table.primes) {
            const u64 q = static_cast<u64>(prime);
            if (q * q > rem) {
                break;
            }
            if (rem % q != 0ULL) {
                continue;
            }
            int count = 0;
            do {
                rem /= q;
                ++count;
            } while (rem % q == 0ULL);
            factors.push_back({q, count});
        }
        if (rem > 1ULL) {
            factors.push_back({rem, 1});
        }

        divisors.clear();
        build_divisors_leq(factors, 0U, 1ULL, root, divisors);

        for (const u64 d : divisors) {
            const u64 e = n / d;
            const u128 pu = static_cast<u128>(p);
            const u128 value = pu * (pu + static_cast<u128>(d)) * (pu + static_cast<u128>(e));

            if (heap.size() < static_cast<std::size_t>(target)) {
                if (active.insert(value).second) {
                    heap.push_back(value);
                    std::push_heap(heap.begin(), heap.end());
                }
                continue;
            }

            const u128 current_max = heap.front();
            if (value >= current_max) {
                continue;
            }
            if (active.find(value) != active.end()) {
                continue;
            }

            std::pop_heap(heap.begin(), heap.end());
            const u128 removed = heap.back();
            heap.pop_back();
            active.erase(removed);

            heap.push_back(value);
            std::push_heap(heap.begin(), heap.end());
            active.insert(value);
        }

        if (heap.size() == static_cast<std::size_t>(target)) {
            const u64 next_p = p + 1ULL;
            const u128 lower_bound = static_cast<u128>(next_p) *
                                     static_cast<u128>(next_p + 1ULL) *
                                     static_cast<u128>(next_p + 1ULL);
            if (lower_bound > heap.front()) {
                return static_cast<u64>(heap.front());
            }
        }
    }
}

bool run_checkpoints() {
    if (nth_alexandrian(6) != 630ULL) {
        std::cerr << "Checkpoint failed for 6th Alexandrian integer" << '\n';
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
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::cout << nth_alexandrian(options.target) << '\n';
    return 0;
}
