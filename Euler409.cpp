#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

constexpr int MOD = 1000000007;

struct Options {
    int n = 10000000;
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

int mod_pow(i64 base, i64 exp) {
    i64 result = 1;
    i64 cur = base % MOD;
    i64 e = exp;
    while (e > 0) {
        if (e & 1LL) {
            result = result * cur % MOD;
        }
        cur = cur * cur % MOD;
        e >>= 1LL;
    }
    return static_cast<int>(result);
}

int W_mod(const int n) {
    const int pow2n = mod_pow(2, n);
    const int total_base = (pow2n - 1 + MOD) % MOD;  // 2^n - 1

    // Total ordered distinct non-zero tuples: T = (2^n-1)_n.
    int T = 1;
    for (int i = 0; i < n; ++i) {
        T = static_cast<int>(1LL * T * ((total_base - i + MOD) % MOD) % MOD);
    }

    const int m = n / 2;

    int fact_n = 1;
    int fact_m = 1;
    for (int i = 1; i <= n; ++i) {
        fact_n = static_cast<int>(1LL * fact_n * i % MOD);
        if (i <= m) {
            fact_m = static_cast<int>(1LL * fact_m * i % MOD);
        }
    }
    const int inv_fact_m = mod_pow(fact_m, MOD - 2);

    const int q_minus_1 = (mod_pow(2, n - 1) - 1 + MOD) % MOD;  // 2^(n-1)-1
    int falling_q = 1;
    for (int i = 0; i < m; ++i) {
        falling_q = static_cast<int>(1LL * falling_q * ((q_minus_1 - i + MOD) % MOD) % MOD);
    }
    int comb = static_cast<int>(1LL * falling_q * inv_fact_m % MOD);  // C(2^(n-1)-1, floor(n/2))

    // Coefficient [x^n] of (1-x)^(2^(n-1)) (1+x)^(2^(n-1)-1):
    // if n even:  (+/-) C(2^(n-1)-1, n/2) with sign (-1)^(n/2)
    // if n odd :  (+/-) C(2^(n-1)-1, (n-1)/2) with sign -(-1)^((n-1)/2)
    if ((n & 1) == 0) {
        if ((m & 1) == 1) {
            comb = (MOD - comb) % MOD;
        }
    } else {
        if (((m) & 1) == 0) {
            comb = (MOD - comb) % MOD;
        }
    }

    const int S = static_cast<int>(1LL * fact_n * comb % MOD);

    const int inv_pow2n = mod_pow(pow2n, MOD - 2);
    const int losing = static_cast<int>(
        1LL * ((T + 1LL * (pow2n - 1 + MOD) % MOD * S) % MOD) * inv_pow2n % MOD);
    const int winning = (T - losing + MOD) % MOD;
    return winning;
}

bool run_checkpoints() {
    if (W_mod(1) != 1) {
        std::cerr << "Checkpoint failed: W(1)\n";
        return false;
    }
    if (W_mod(2) != 6) {
        std::cerr << "Checkpoint failed: W(2)\n";
        return false;
    }
    if (W_mod(3) != 168) {
        std::cerr << "Checkpoint failed: W(3)\n";
        return false;
    }
    if (W_mod(5) != 19764360) {
        std::cerr << "Checkpoint failed: W(5)\n";
        return false;
    }
    if (W_mod(100) != 384777056) {
        std::cerr << "Checkpoint failed: W(100)\n";
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

    std::cout << W_mod(options.n) << '\n';
    return 0;
}
