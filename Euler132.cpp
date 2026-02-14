#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int wanted = 40;
    u64 exponent = 1000000000ULL;
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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

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
        if (parse_int_after_prefix(arg, "--wanted=", options.wanted) ||
            parse_u64_after_prefix(arg, "--exponent=", options.exponent)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.wanted >= 1;
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 repunit_mod(u64 n, u64 mod) {
    if (mod == 1) {
        return 0;
    }

    u64 pow10 = 1 % mod;
    u64 repunit = 0;
    const u64 ten_mod = 10 % mod;

    int msb = 63;
    while (msb > 0 && ((n >> msb) & 1ULL) == 0ULL) {
        --msb;
    }

    for (int bit = msb; bit >= 0; --bit) {
        repunit = (mul_mod(repunit, pow10, mod) + repunit) % mod;
        pow10 = mul_mod(pow10, pow10, mod);

        if (((n >> bit) & 1ULL) != 0ULL) {
            repunit = (mul_mod(repunit, ten_mod, mod) + 1ULL) % mod;
            pow10 = mul_mod(pow10, ten_mod, mod);
        }
    }

    return repunit;
}

bool is_prime(int n, const std::vector<int>& primes) {
    if (n <= 1) {
        return false;
    }
    for (int p : primes) {
        if (static_cast<std::int64_t>(p) * p > n) {
            break;
        }
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

bool divides_repunit(int p, u64 exponent) {
    return repunit_mod(exponent, static_cast<u64>(p)) == 0ULL;
}

u64 solve(const int wanted, const u64 exponent) {
    int found = 0;
    u64 sum = 0;
    std::vector<int> primes;
    primes.reserve(4096);
    primes.push_back(2);

    if (wanted <= 0) {
        return 0;
    }

    for (int p = 3; found < wanted; p += 2) {
        if (!is_prime(p, primes)) {
            continue;
        }
        primes.push_back(p);

        if (p == 5) {
            continue;
        }
        if (divides_repunit(p, exponent)) {
            ++found;
            sum += static_cast<u64>(p);
        }
    }

    return sum;
}

bool run_checkpoints() {
    if (solve(4, 10ULL) != 9414ULL) {
        std::cerr << "Checkpoint failed for first 4 factors of R(10)" << '\n';
        return false;
    }
    if (solve(4, 6ULL) != 34ULL) {
        std::cerr << "Checkpoint failed for first 4 factors of R(10^6)" << '\n';
        return false;
    }
    if (divides_repunit(3, 10ULL) || !divides_repunit(3, 3ULL)) {
        std::cerr << "Checkpoint failed for prime 3 divisibility behavior" << '\n';
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

    std::cout << solve(options.wanted, options.exponent) << '\n';
    return 0;
}
