#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int p = 1009;
    int q = 3643;
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
            parse_int_after_prefix(arg, "--q=", options.q)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.p > 2 && options.q > 2 && options.p != options.q;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    base %= mod;
    u64 result = 1 % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

i64 unconcealed_count_formula(const int e, const int p, const int q) {
    const i64 gp = std::gcd(e - 1, p - 1);
    const i64 gq = std::gcd(e - 1, q - 1);
    return (gp + 1) * (gq + 1);
}

i64 unconcealed_count_bruteforce(const int e, const int p, const int q) {
    const int n = p * q;
    i64 count = 0;
    for (int m = 0; m < n; ++m) {
        if (pow_mod(static_cast<u64>(m), static_cast<u64>(e), static_cast<u64>(n)) == static_cast<u64>(m)) {
            ++count;
        }
    }
    return count;
}

i64 solve(const int p, const int q) {
    const int phi = (p - 1) * (q - 1);
    i64 best = -1;
    i64 sum = 0;

    for (int e = 2; e < phi; ++e) {
        if (std::gcd(e, phi) != 1) {
            continue;
        }

        const i64 cnt = unconcealed_count_formula(e, p, q);
        if (best < 0 || cnt < best) {
            best = cnt;
            sum = e;
        } else if (cnt == best) {
            sum += e;
        }
    }

    return sum;
}

i64 brute_solve(const int p, const int q) {
    const int phi = (p - 1) * (q - 1);
    i64 best = -1;
    i64 sum = 0;

    for (int e = 2; e < phi; ++e) {
        if (std::gcd(e, phi) != 1) {
            continue;
        }

        const i64 cnt = unconcealed_count_bruteforce(e, p, q);
        if (best < 0 || cnt < best) {
            best = cnt;
            sum = e;
        } else if (cnt == best) {
            sum += e;
        }
    }

    return sum;
}

bool run_checkpoints() {
    if (unconcealed_count_formula(181, 19, 37) != 703LL) {
        std::cerr << "Checkpoint failed for stated sample e=181, p=19, q=37" << '\n';
        return false;
    }

    for (int e = 2; e < (5 - 1) * (11 - 1); ++e) {
        if (std::gcd(e, (5 - 1) * (11 - 1)) != 1) {
            continue;
        }
        if (unconcealed_count_formula(e, 5, 11) != unconcealed_count_bruteforce(e, 5, 11)) {
            std::cerr << "Checkpoint failed for formula-vs-brute at p=5, q=11" << '\n';
            return false;
        }
    }

    if (solve(19, 37) != brute_solve(19, 37)) {
        std::cerr << "Checkpoint failed for full brute cross-check at p=19, q=37" << '\n';
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

    std::cout << solve(options.p, options.q) << '\n';
    return 0;
}
