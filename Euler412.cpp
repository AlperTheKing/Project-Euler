#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

constexpr int MOD = 76543217;

struct Options {
    int m = 10000;
    int n = 5000;
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
        if (parse_int_after_prefix(arg, "--m=", options.m) ||
            parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.m >= 1 && options.n >= 0 && options.n < options.m;
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

int LC_mod(const int m, const int n) {
    const int p = m - n;
    const i64 cells = 1LL * m * m - 1LL * n * n;

    const int max_small = 2 * m - 1;
    std::vector<int> fact_small(static_cast<std::size_t>(max_small + 1), 1);
    for (int i = 1; i <= max_small; ++i) {
        fact_small[static_cast<std::size_t>(i)] =
            static_cast<int>(1LL * fact_small[static_cast<std::size_t>(i - 1)] * i % MOD);
    }

    int fact_cells = 1;
    for (i64 i = 1; i <= cells; ++i) {
        fact_cells = static_cast<int>(1LL * fact_cells * (i % MOD) % MOD);
    }

    // Hook-product decomposition for shape λ = (m repeated p times, p repeated n times):
    // P = P_A * P_B^2, where
    // P_A = prod_{t=m+n}^{2m-1} t! / prod_{s=2n}^{m+n-1} s!
    // P_B = prod_{t=n}^{m-1} t! / prod_{u=0}^{p-1} u!
    int pa_num = 1;
    for (int t = m + n; t <= 2 * m - 1; ++t) {
        pa_num = static_cast<int>(1LL * pa_num * fact_small[static_cast<std::size_t>(t)] % MOD);
    }
    int pa_den = 1;
    for (int s = 2 * n; s <= m + n - 1; ++s) {
        pa_den = static_cast<int>(1LL * pa_den * fact_small[static_cast<std::size_t>(s)] % MOD);
    }
    const int pa = static_cast<int>(1LL * pa_num * mod_pow(pa_den, MOD - 2) % MOD);

    int pb_num = 1;
    for (int t = n; t <= m - 1; ++t) {
        pb_num = static_cast<int>(1LL * pb_num * fact_small[static_cast<std::size_t>(t)] % MOD);
    }
    int pb_den = 1;
    for (int u = 0; u <= p - 1; ++u) {
        pb_den = static_cast<int>(1LL * pb_den * fact_small[static_cast<std::size_t>(u)] % MOD);
    }
    const int pb = static_cast<int>(1LL * pb_num * mod_pow(pb_den, MOD - 2) % MOD);

    int hooks = pa;
    hooks = static_cast<int>(1LL * hooks * pb % MOD);
    hooks = static_cast<int>(1LL * hooks * pb % MOD);

    return static_cast<int>(1LL * fact_cells * mod_pow(hooks, MOD - 2) % MOD);
}

cpp_int LC_exact_small(const int m, const int n) {
    const int p = m - n;
    const int rows = m;
    std::vector<int> row_len(static_cast<std::size_t>(rows), p);
    for (int i = 0; i < p; ++i) {
        row_len[static_cast<std::size_t>(i)] = m;
    }

    std::vector<int> col_len(static_cast<std::size_t>(m), p);
    for (int j = 0; j < p; ++j) {
        col_len[static_cast<std::size_t>(j)] = m;
    }

    int cells = 0;
    for (int i = 0; i < rows; ++i) {
        cells += row_len[static_cast<std::size_t>(i)];
    }

    cpp_int numer = 1;
    for (int i = 1; i <= cells; ++i) {
        numer *= i;
    }

    cpp_int denom = 1;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < row_len[static_cast<std::size_t>(i)]; ++j) {
            const int hook =
                (row_len[static_cast<std::size_t>(i)] - j) + (col_len[static_cast<std::size_t>(j)] - i) - 1;
            denom *= hook;
        }
    }

    return numer / denom;
}

bool run_checkpoints() {
    if (LC_exact_small(3, 0) != 42) {
        std::cerr << "Checkpoint failed: LC(3,0)\n";
        return false;
    }
    if (LC_exact_small(5, 3) != cpp_int("250250")) {
        std::cerr << "Checkpoint failed: LC(5,3)\n";
        return false;
    }
    if (LC_exact_small(6, 3) != cpp_int("406029023400")) {
        std::cerr << "Checkpoint failed: LC(6,3)\n";
        return false;
    }
    if (LC_mod(10, 5) != 61251715) {
        std::cerr << "Checkpoint failed: LC(10,5) mod 76543217\n";
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

    std::cout << LC_mod(options.m, options.n) << '\n';
    return 0;
}
