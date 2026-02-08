#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;

struct Options {
    i64 reflections = 12017639147LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--reflections=", options.reflections)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.reflections >= 1;
}

std::vector<i64> distinct_prime_factors(i64 n) {
    std::vector<i64> factors;
    for (i64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        factors.push_back(p);
        while (n % p == 0) {
            n /= p;
        }
    }
    if (n > 1) {
        factors.push_back(n);
    }
    return factors;
}

i64 count_residue_upto(const i64 m, const int residue) {
    if (m <= 0) {
        return 0;
    }

    i64 first = residue;
    if (first == 0) {
        first = 3;
    }
    if (first > m) {
        return 0;
    }

    return 1 + (m - first) / 3;
}

i64 count_multiples_with_mod2(const i64 n, const i64 d) {
    if (d % 3 == 0) {
        return 0;
    }
    const i64 m = (n - 1) / d;
    const int inv_mod3 = (d % 3 == 1) ? 1 : 2;
    const int residue_t = (2 * inv_mod3) % 3;
    return count_residue_upto(m, residue_t);
}

i64 solve(const i64 reflections) {
    const i64 n = (reflections + 3) / 2;
    const std::vector<i64> factors = distinct_prime_factors(n);

    i64 answer = 0;
    const auto dfs = [&](const auto& self, const std::size_t index, const i64 product, const int parity) -> void {
        if (index == factors.size()) {
            const i64 cnt = count_multiples_with_mod2(n, product);
            if ((parity & 1) == 0) {
                answer += cnt;
            } else {
                answer -= cnt;
            }
            return;
        }

        self(self, index + 1, product, parity);

        const i64 p = factors[index];
        if (product > (n - 1) / p) {
            return;
        }
        self(self, index + 1, product * p, parity + 1);
    };

    dfs(dfs, 0U, 1, 0);
    return answer;
}

bool run_checkpoints() {
    if (solve(11) != 2) {
        std::cerr << "Checkpoint failed for 11 reflections" << '\n';
        return false;
    }
    if (solve(1000001) != 80840) {
        std::cerr << "Checkpoint failed for 1000001 reflections" << '\n';
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

    std::cout << solve(options.reflections) << '\n';
    return 0;
}
