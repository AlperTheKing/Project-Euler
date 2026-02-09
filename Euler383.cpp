#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using i64 = long long;

struct Options {
    i64 n = 1000000000000000000LL;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(ch - '0');
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

int carry_out(const int digit, const int carry_in) {
    return (2 * digit + carry_in) / 5;
}

struct Key {
    i64 x;
    int k;
    int carry_in;

    bool operator==(const Key& other) const {
        return x == other.x && k == other.k && carry_in == other.carry_in;
    }
};

struct KeyHash {
    std::size_t operator()(const Key& key) const {
        std::size_t h = std::hash<i64>{}(key.x);
        h ^= std::hash<int>{}(key.k) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(key.carry_in) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

class Solver {
   public:
    i64 solve(const i64 n) {
        memo_.clear();
        i64 answer = 0;
        i64 power_of_five = 5;

        for (int a = 1; power_of_five <= n; ++a) {
            answer += count_non_multiples_of_five(n / power_of_five, a - 1);
            if (power_of_five > n / 5) {
                break;
            }
            power_of_five *= 5;
        }

        return answer;
    }

   private:
    std::unordered_map<Key, i64, KeyHash> memo_;

    i64 count_with_initial_carry(const i64 x, const int k, const int carry_in) {
        if (x < 0 || k < 0) {
            return 0;
        }
        if (x == 0) {
            // Only n=0 is present, and it contributes zero further carries.
            return 1;
        }

        const Key key{x, k, carry_in};
        const auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        i64 ways = 0;
        for (int digit = 0; digit <= 4; ++digit) {
            if (x < digit) {
                break;
            }
            const i64 q = (x - digit) / 5;
            const int next_carry = carry_out(digit, carry_in);
            ways += count_with_initial_carry(q, k - next_carry, next_carry);
        }

        memo_[key] = ways;
        return ways;
    }

    i64 count_non_multiples_of_five(const i64 x, const int k) {
        if (x <= 0 || k < 0) {
            return 0;
        }
        i64 ways = 0;
        for (int digit = 1; digit <= 4; ++digit) {
            if (x < digit) {
                break;
            }
            const i64 q = (x - digit) / 5;
            const int next_carry = carry_out(digit, 0);
            ways += count_with_initial_carry(q, k - next_carry, next_carry);
        }
        return ways;
    }
};

int v5_factorial(i64 n) {
    int total = 0;
    while (n > 0) {
        n /= 5;
        total += static_cast<int>(n);
    }
    return total;
}

i64 brute_force_count(const i64 n) {
    i64 count = 0;
    for (i64 i = 1; i <= n; ++i) {
        const int lhs = v5_factorial(2 * i - 1);
        const int rhs = 2 * v5_factorial(i);
        if (lhs < rhs) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    Solver solver;

    if (solver.solve(1000) != 68LL) {
        std::cerr << "Checkpoint failed: T5(10^3)\n";
        return false;
    }

    if (solver.solve(1000000000LL) != 2408210LL) {
        std::cerr << "Checkpoint failed: T5(10^9)\n";
        return false;
    }

    const i64 brute_n = 20000;
    const i64 expected = brute_force_count(brute_n);
    const i64 got = solver.solve(brute_n);
    if (got != expected) {
        std::cerr << "Checkpoint failed against brute force at n=" << brute_n
                  << " (expected " << expected << ", got " << got << ")\n";
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

    Solver solver;
    std::cout << solver.solve(options.n) << '\n';
    return 0;
}
