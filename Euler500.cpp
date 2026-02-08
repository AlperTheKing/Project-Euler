#include <cmath>
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int power = 500'500;
    u64 mod = 500'500'507ULL;
    bool run_checkpoints = true;
};

struct HeapNode {
    long double log_value;
    u64 mod_value;
};

struct HeapCmp {
    bool operator()(const HeapNode& lhs, const HeapNode& rhs) const {
        return lhs.log_value > rhs.log_value;
    }
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--power=", options.power) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.power >= 1 && options.mod >= 2ULL;
}

int nth_prime_upper_bound(const int n) {
    if (n < 6) {
        return 20;
    }
    const long double x = static_cast<long double>(n);
    const long double estimate = x * (std::log(x) + std::log(std::log(x))) + 10.0L;
    return static_cast<int>(estimate) + 100;
}

std::vector<int> first_n_primes(const int n) {
    int limit = nth_prime_upper_bound(n);
    while (true) {
        std::vector<bool> is_composite(static_cast<std::size_t>(limit + 1), false);
        std::vector<int> primes;
        primes.reserve(static_cast<std::size_t>(n + 10));
        for (int i = 2; i <= limit; ++i) {
            if (!is_composite[static_cast<std::size_t>(i)]) {
                primes.push_back(i);
                if (i <= limit / i) {
                    for (int j = i * i; j <= limit; j += i) {
                        is_composite[static_cast<std::size_t>(j)] = true;
                    }
                }
            }
        }
        if (static_cast<int>(primes.size()) >= n) {
            primes.resize(static_cast<std::size_t>(n));
            return primes;
        }
        limit *= 2;
    }
}

u64 solve(const int power, const u64 mod) {
    const std::vector<int> primes = first_n_primes(power);
    std::priority_queue<HeapNode, std::vector<HeapNode>, HeapCmp> heap;
    heap = {};
    for (int p : primes) {
        heap.push(HeapNode{std::log(static_cast<long double>(p)), static_cast<u64>(p) % mod});
    }

    u64 ans = 1ULL;
    for (int step = 0; step < power; ++step) {
        HeapNode cur = heap.top();
        heap.pop();
        ans = static_cast<u64>((static_cast<u128>(ans) * cur.mod_value) % mod);
        const u64 squared_mod = static_cast<u64>((static_cast<u128>(cur.mod_value) * cur.mod_value) % mod);
        heap.push(HeapNode{cur.log_value * 2.0L, squared_mod});
    }
    return ans;
}

u64 brute_smallest_with_two_power_divisors(const int power) {
    const u64 target = 1ULL << power;
    for (u64 n = 1ULL;; ++n) {
        u64 x = n;
        u64 dcount = 1ULL;
        for (u64 p = 2ULL; p * p <= x; ++p) {
            if (x % p != 0ULL) {
                continue;
            }
            int e = 0;
            while (x % p == 0ULL) {
                x /= p;
                ++e;
            }
            dcount *= static_cast<u64>(e + 1);
        }
        if (x > 1ULL) {
            dcount *= 2ULL;
        }
        if (dcount == target) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (solve(4, 1'000'000'007ULL) != 120ULL) {
        std::cerr << "Checkpoint failed: smallest with 16 divisors is 120" << '\n';
        return false;
    }
    if (solve(8, 1'000'000'007ULL) != brute_smallest_with_two_power_divisors(8)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for power=8" << '\n';
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
    std::cout << solve(options.power, options.mod) << '\n';
    return 0;
}
