#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 target_num = 15499;
    u64 target_den = 94744;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--target-num=", options.target_num) ||
            parse_u64_after_prefix(arg, "--target-den=", options.target_den)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.target_den != 0;
}

bool resilience_less_than(const u64 n, const u64 phi_n, const u64 target_num, const u64 target_den) {
    return static_cast<u128>(phi_n) * static_cast<u128>(target_den) <
           static_cast<u128>(target_num) * static_cast<u128>(n - 1);
}

u64 solve(const u64 target_num, const u64 target_den) {
    const std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    u64 best = ~0ULL;

    const auto dfs = [&](const auto& self, const std::size_t prime_index, const int max_exp,
                         const u64 n, const u64 phi_n) -> void {
        if (n >= best) {
            return;
        }

        if (n > 1 && resilience_less_than(n, phi_n, target_num, target_den)) {
            best = n;
            return;
        }

        if (prime_index >= primes.size()) {
            return;
        }

        const u64 p = static_cast<u64>(primes[prime_index]);
        u64 n_mul = n;
        u64 phi_mul = phi_n;

        for (int exp = 1; exp <= max_exp; ++exp) {
            if (n_mul > best / p) {
                break;
            }
            n_mul *= p;
            if (exp == 1) {
                phi_mul *= (p - 1);
            } else {
                phi_mul *= p;
            }
            self(self, prime_index + 1, exp, n_mul, phi_mul);
        }
    };

    dfs(dfs, 0U, 64, 1, 1);
    return best;
}

bool run_checkpoints() {
    if (solve(4, 10) != 12ULL) {
        std::cerr << "Checkpoint failed for threshold 4/10" << '\n';
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

    std::cout << solve(options.target_num, options.target_den) << '\n';
    return 0;
}
