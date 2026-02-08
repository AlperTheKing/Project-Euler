#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int p = 61;
    int q = 10000000;
    int exponent = 10;
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
        if (parse_int_after_prefix(arg, "--p=", options.p) ||
            parse_int_after_prefix(arg, "--q=", options.q) ||
            parse_int_after_prefix(arg, "--exponent=", options.exponent)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.p >= 2 && options.q >= 0 && options.exponent >= 1 && options.exponent <= 24;
}

u64 pow_u64(const u64 base, const int exp) {
    u64 value = 1;
    for (int i = 0; i < exp; ++i) {
        value *= base;
    }
    return value;
}

u64 solve(const int p, const int q, const int exponent) {
    const u64 mod = pow_u64(static_cast<u64>(p), exponent);
    const int first_len = std::min(exponent, q + 1);
    std::vector<u64> prefix(static_cast<std::size_t>(first_len + 1), 0ULL);

    u64 s = 290797ULL;
    u64 total = 0ULL;
    for (int n = 0; n <= q; ++n) {
        const u64 t = s % static_cast<u64>(p);
        total += t;
        if (n < first_len) {
            prefix[static_cast<std::size_t>(n + 1)] = prefix[static_cast<std::size_t>(n)] + t;
        }
        s = static_cast<u64>((static_cast<u128>(s) * s) % 50515093ULL);
    }

    u64 answer = 0ULL;
    u64 p_pow = 1ULL;
    for (int j = 0; j < exponent; ++j) {
        u64 pref = 0ULL;
        if (j + 1 <= first_len) {
            pref = prefix[static_cast<std::size_t>(j + 1)];
        } else if (first_len == q + 1) {
            pref = total;
        } else {
            pref = prefix[static_cast<std::size_t>(first_len)];
        }
        const u64 suffix = total - pref;
        answer = static_cast<u64>((static_cast<u128>(answer) +
                                   static_cast<u128>(p_pow) * (suffix % mod)) %
                                  mod);
        p_pow = static_cast<u64>((static_cast<u128>(p_pow) * static_cast<u64>(p)) % mod);
    }
    return answer;
}

u64 brute_small(const int p, const int q, const int exponent) {
    const u64 mod = pow_u64(static_cast<u64>(p), exponent);
    u64 s = 290797ULL;
    cpp_int n_value = 0;
    cpp_int p_pow = 1;
    for (int i = 0; i <= q; ++i) {
        const int t = static_cast<int>(s % static_cast<u64>(p));
        n_value += cpp_int(t) * p_pow;
        p_pow *= p;
        s = static_cast<u64>((static_cast<u128>(s) * s) % 50515093ULL);
    }

    cpp_int valuation = 0;
    cpp_int current = n_value;
    while (current > 0) {
        current /= p;
        valuation += current;
    }
    return static_cast<u64>(valuation % mod);
}

bool run_checkpoints() {
    if (solve(3, 10000, 20) != 624955285ULL) {
        std::cerr << "Checkpoint failed for stated sample NF(3,10000) mod 3^20" << '\n';
        return false;
    }
    if (solve(5, 200, 8) != brute_small(5, 200, 8)) {
        std::cerr << "Checkpoint failed for brute cross-check at p=5 q=200" << '\n';
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
    std::cout << solve(options.p, options.q, options.exponent) << '\n';
    return 0;
}
