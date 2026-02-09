#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

constexpr u64 MOD = 1000000007ULL;

struct Options {
    int limit = 50000000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    u64 cur = base % MOD;
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = static_cast<u64>((__uint128_t)result * cur % MOD);
        }
        cur = static_cast<u64>((__uint128_t)cur * cur % MOD);
        e >>= 1ULL;
    }
    return result;
}

std::vector<bool> sieve_primes(int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; 1LL * p * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= n; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }
    return is_prime;
}

std::vector<u64> modular_inverses(int n) {
    std::vector<u64> inv(static_cast<std::size_t>(n + 1), 0ULL);
    if (n >= 1) {
        inv[1] = 1ULL;
    }
    for (int i = 2; i <= n; ++i) {
        inv[static_cast<std::size_t>(i)] =
            (MOD - (MOD / static_cast<u64>(i)) * inv[static_cast<std::size_t>(MOD % static_cast<u64>(i))] % MOD) % MOD;
    }
    return inv;
}

u64 solve_mod(const int limit) {
    const std::vector<bool> is_prime = sieve_primes(limit + 1);
    const std::vector<u64> inv = modular_inverses(limit + 2);
    const u64 inv5 = mod_pow(5ULL, MOD - 2ULL);

    int t = 0;        // t = pi(n)
    u64 p = 1ULL;     // sum_{k=0..t} B_n(k)
    u64 bt = 1ULL;    // B_n(t)
    u64 bt1 = 0ULL;   // B_n(t+1)

    u64 answer = 0ULL;

    for (int n = 1; n <= limit; ++n) {
        const u64 c_n = 6ULL * p % MOD;
        answer += c_n;
        if (answer >= MOD) {
            answer -= MOD;
        }

        if (n == limit) {
            break;
        }

        const bool next_is_prime = is_prime[static_cast<std::size_t>(n + 1)];

        u64 p_next = 0ULL;
        u64 bt_next = 0ULL;
        u64 bt1_next = 0ULL;

        if (!next_is_prime) {
            p_next = (6ULL * p + MOD - bt) % MOD;

            u64 bt_minus_1 = 0ULL;
            if (t > 0) {
                bt_minus_1 = bt * (5ULL * static_cast<u64>(t) % MOD) % MOD;
                bt_minus_1 = bt_minus_1 * inv[static_cast<std::size_t>(n - t)] % MOD;
            }

            bt_next = (5ULL * bt + bt_minus_1) % MOD;
            bt1_next = (5ULL * bt1 + bt) % MOD;
        } else {
            p_next = (6ULL * p + 5ULL * bt1) % MOD;

            u64 bt2 = 0ULL;  // B_n(t+2)
            if (t + 2 <= n - 1) {
                bt2 = bt1 * static_cast<u64>(n - t - 2) % MOD;
                bt2 = bt2 * inv5 % MOD;
                bt2 = bt2 * inv[static_cast<std::size_t>(t + 2)] % MOD;
            }

            bt_next = (5ULL * bt1 + bt) % MOD;
            bt1_next = (5ULL * bt2 + bt1) % MOD;
            ++t;
        }

        p = p_next;
        bt = bt_next;
        bt1 = bt1_next;
    }

    return answer;
}

cpp_int comb_small(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    const int r = std::min(k, n - k);
    cpp_int ans = 1;
    for (int i = 1; i <= r; ++i) {
        ans *= (n - r + i);
        ans /= i;
    }
    return ans;
}

cpp_int pow_cpp_int(cpp_int base, int exp) {
    cpp_int result = 1;
    int e = exp;
    while (e > 0) {
        if (e & 1) {
            result *= base;
        }
        base *= base;
        e >>= 1;
    }
    return result;
}

cpp_int c_exact(int n) {
    const std::vector<bool> is_prime = sieve_primes(n);
    int pi_n = 0;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            ++pi_n;
        }
    }

    cpp_int sum = 0;
    for (int k = 0; k <= pi_n; ++k) {
        sum += comb_small(n - 1, k) * pow_cpp_int(cpp_int(5), n - 1 - k);
    }
    return sum * 6;
}

bool run_checkpoints() {
    if (c_exact(3) != 216) {
        std::cerr << "Checkpoint failed: C(3)\n";
        return false;
    }
    if (c_exact(4) != 1290) {
        std::cerr << "Checkpoint failed: C(4)\n";
        return false;
    }
    if (c_exact(11) != cpp_int("361912500")) {
        std::cerr << "Checkpoint failed: C(11)\n";
        return false;
    }
    if (c_exact(24) != cpp_int("4727547363281250000")) {
        std::cerr << "Checkpoint failed: C(24)\n";
        return false;
    }
    if (solve_mod(50) != 832833871ULL) {
        std::cerr << "Checkpoint failed: S(50) mod 1e9+7\n";
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

    std::cout << solve_mod(options.limit) << '\n';
    return 0;
}
