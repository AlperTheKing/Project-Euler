#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int set_size = 5;
    int prime_limit = 10000;
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
        if (parse_int_after_prefix(arg, "--set-size=", options.set_size)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--prime-limit=", options.prime_limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.set_size >= 2 && options.prime_limit >= 100;
}

u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<__uint128_t>(a) * b) % mod);
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1ULL;
    base %= mod;

    while (exp > 0) {
        if (exp & 1ULL) {
            result = mod_mul(result, base, mod);
        }
        base = mod_mul(base, base, mod);
        exp >>= 1ULL;
    }

    return result;
}

bool is_prime_u64(u64 n) {
    if (n < 2ULL) {
        return false;
    }
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0ULL) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    // Deterministic bases for 64-bit integers.
    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0ULL) {
            continue;
        }

        u64 x = mod_pow(a, d, n);
        if (x == 1ULL || x == n - 1ULL) {
            continue;
        }

        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mod_mul(x, x, n);
            if (x == n - 1ULL) {
                witness = false;
                break;
            }
        }

        if (witness) {
            return false;
        }
    }

    return true;
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    is_prime[0] = false;
    is_prime[1] = false;

    for (int p = 2; p <= limit / p; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] && p != 2 && p != 5) {
            primes.push_back(p);
        }
    }

    return primes;
}

u64 concat_numbers(const int a, const int b) {
    u64 pow10 = 10ULL;
    int t = b;
    while (t >= 10) {
        t /= 10;
        pow10 *= 10ULL;
    }
    return static_cast<u64>(a) * pow10 + static_cast<u64>(b);
}

int solve(const int set_size, const int prime_limit) {
    const std::vector<int> primes = sieve_primes(prime_limit);
    const int n = static_cast<int>(primes.size());

    std::vector<std::vector<unsigned char>> compatible(
        static_cast<std::size_t>(n), std::vector<unsigned char>(static_cast<std::size_t>(n), 0));

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const u64 ab = concat_numbers(primes[i], primes[j]);
            const u64 ba = concat_numbers(primes[j], primes[i]);
            if (is_prime_u64(ab) && is_prime_u64(ba)) {
                compatible[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 1;
                compatible[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = 1;
            }
        }
    }

    int best_sum = std::numeric_limits<int>::max();

    std::vector<int> chosen;
    std::vector<int> all_indices(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        all_indices[static_cast<std::size_t>(i)] = i;
    }

    auto dfs = [&](auto&& self, std::vector<int>& candidates, int current_sum) -> void {
        if (static_cast<int>(chosen.size()) == set_size) {
            if (current_sum < best_sum) {
                best_sum = current_sum;
            }
            return;
        }

        const int need = set_size - static_cast<int>(chosen.size());
        if (static_cast<int>(candidates.size()) < need) {
            return;
        }

        for (std::size_t idx = 0; idx < candidates.size(); ++idx) {
            const int p_index = candidates[idx];
            const int next_sum = current_sum + primes[static_cast<std::size_t>(p_index)];
            if (next_sum >= best_sum) {
                continue;
            }

            chosen.push_back(p_index);

            std::vector<int> next_candidates;
            next_candidates.reserve(candidates.size());

            for (std::size_t j = idx + 1; j < candidates.size(); ++j) {
                const int q_index = candidates[j];
                bool ok = true;
                for (const int c_index : chosen) {
                    if (c_index == q_index) {
                        continue;
                    }
                    if (!compatible[static_cast<std::size_t>(c_index)][static_cast<std::size_t>(q_index)]) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    next_candidates.push_back(q_index);
                }
            }

            self(self, next_candidates, next_sum);
            chosen.pop_back();
        }
    };

    dfs(dfs, all_indices, 0);
    return best_sum;
}

bool run_checkpoints() {
    if (solve(4, 1000) != 792) {
        std::cerr << "Checkpoint failed for set_size=4" << '\n';
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

    std::cout << solve(options.set_size, options.prime_limit) << '\n';
    return 0;
}
