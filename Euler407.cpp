#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = long long;

struct Options {
    int limit = 10000000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

std::vector<int> build_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (const int p : primes) {
            const i64 x = 1LL * i * p;
            if (x > n) {
                break;
            }
            spf[static_cast<std::size_t>(x)] = p;
            if (p == spf[static_cast<std::size_t>(i)]) {
                break;
            }
        }
    }
    return spf;
}

int inv_mod(int a, int mod) {
    int b = mod;
    int u = 1;
    int v = 0;
    while (b != 0) {
        const int t = a / b;
        a -= t * b;
        std::swap(a, b);
        u -= t * v;
        std::swap(u, v);
    }
    if (u < 0) {
        u += mod;
    }
    return u;
}

u64 sum_M(const int limit) {
    const std::vector<int> spf = build_spf(limit);

    std::vector<int> prime_powers(10, 0);
    std::vector<int> subset_prod(1 << 10, 1);

    u64 ans = 0ULL;
    for (int n = 1; n <= limit; ++n) {
        if (n == 1) {
            continue;
        }

        int x = n;
        int k = 0;
        while (x > 1) {
            const int p = spf[static_cast<std::size_t>(x)];
            int pe = 1;
            while (x % p == 0) {
                x /= p;
                pe *= p;
            }
            prime_powers[static_cast<std::size_t>(k++)] = pe;
        }

        if (k == 1) {
            ans += 1ULL;
            continue;
        }

        const int full = (1 << k);
        subset_prod[0] = 1;
        for (int mask = 1; mask < full; ++mask) {
            const int lsb = mask & -mask;
            const int bit = __builtin_ctz(static_cast<unsigned>(lsb));
            subset_prod[mask] = subset_prod[mask ^ lsb] * prime_powers[static_cast<std::size_t>(bit)];
        }

        int best = 1;
        for (int mask = 1; mask < full - 1; ++mask) {
            const int u = subset_prod[mask];
            const int v = n / u;
            const int inv = inv_mod(u % v, v);
            const int a = static_cast<int>((1LL * u * inv) % n);
            if (a > best) {
                best = a;
            }
        }

        ans += static_cast<u64>(best);
    }

    return ans;
}

bool run_checkpoints() {
    // From statement: M(6)=4.
    const u64 s6 = sum_M(6) - sum_M(5);
    if (s6 != 4ULL) {
        std::cerr << "Checkpoint failed: M(6)\n";
        return false;
    }
    // Small consistency: M(1)=0 and M(prime)=1.
    const u64 s2 = sum_M(2);
    if (s2 != 1ULL) {
        std::cerr << "Checkpoint failed: sum up to 2\n";
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

    std::cout << sum_M(options.limit) << '\n';
    return 0;
}
