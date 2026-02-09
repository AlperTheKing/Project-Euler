#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<__int128>(result) * base) % mod);
        }
        base = static_cast<u64>((static_cast<__int128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return result;
}

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2ULL) == 0ULL) {
        return n == 2;
    }
    for (u64 d = 3; d * d <= n; d += 2) {
        if ((n % d) == 0ULL) {
            return false;
        }
    }
    return true;
}

std::vector<u64> prime_factors_distinct(u64 n) {
    std::vector<u64> factors;
    if ((n % 2ULL) == 0ULL) {
        factors.push_back(2ULL);
        while ((n % 2ULL) == 0ULL) {
            n /= 2ULL;
        }
    }
    for (u64 d = 3; d * d <= n; d += 2) {
        if ((n % d) == 0ULL) {
            factors.push_back(d);
            while ((n % d) == 0ULL) {
                n /= d;
            }
        }
    }
    if (n > 1ULL) {
        factors.push_back(n);
    }
    return factors;
}

bool is_full_reptend_prime_base10(const u64 p) {
    if (!is_prime(p) || p == 2ULL || p == 5ULL) {
        return false;
    }
    const u64 phi = p - 1ULL;
    const std::vector<u64> factors = prime_factors_distinct(phi);
    for (const u64 q : factors) {
        if (mod_pow(10ULL, phi / q, p) == 1ULL) {
            return false;
        }
    }
    return true;
}

u64 find_prime_p() {
    // First 11 digits of repetend of 1/p are 00000000137.
    const u64 low = 100000000000ULL / 138ULL + 1ULL;
    const u64 high = 100000000000ULL / 137ULL;

    // Last 5 digits of repetend are 56789.
    // If N is repetend integer, p*N = 10^(p-1)-1 => p*N == -1 (mod 100000).
    // Therefore p*56789 == -1 (mod 100000), i.e. p == 9891 (mod 100000).
    const u64 mod = 100000ULL;
    const u64 residue = 9891ULL;

    u64 start = low;
    const u64 rem = start % mod;
    if (rem <= residue) {
        start += residue - rem;
    } else {
        start += mod - (rem - residue);
    }

    u64 found = 0;
    for (u64 p = start; p <= high; p += mod) {
        if (100000000000ULL / p != 137ULL) {
            continue;
        }
        if (!is_full_reptend_prime_base10(p)) {
            continue;
        }
        if ((p * 56789ULL + 1ULL) % mod != 0ULL) {
            continue;
        }

        if (found != 0ULL) {
            std::cerr << "More than one candidate prime found\n";
            return 0ULL;
        }
        found = p;
    }

    return found;
}

bool run_checkpoints() {
    if (is_full_reptend_prime_base10(7ULL) == false) {
        std::cerr << "Checkpoint failed: p=7 should be full reptend in base 10\n";
        return false;
    }

    const u64 p = find_prime_p();
    if (p == 0ULL) {
        std::cerr << "Checkpoint failed: no unique prime found\n";
        return false;
    }

    if (p != 729809891ULL) {
        std::cerr << "Checkpoint failed: unexpected prime value\n";
        return false;
    }

    if ((9ULL * (p - 1ULL)) / 2ULL != 3284144505ULL) {
        std::cerr << "Checkpoint failed: digit sum formula mismatch\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const u64 p = find_prime_p();
    const u64 answer = (9ULL * (p - 1ULL)) / 2ULL;
    std::cout << answer << '\n';
    return 0;
}
