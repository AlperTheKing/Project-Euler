#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    int k_max = 30;
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
        if (parse_int_after_prefix(arg, "--k-max=", options.k_max)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.k_max >= 1;
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

void radix_sort_u64(std::vector<u64>& values) {
    if (values.size() <= 1U) {
        return;
    }

    std::vector<u64> tmp(values.size());
    constexpr int RADIX_BITS = 16;
    constexpr int RADIX_SIZE = 1 << RADIX_BITS;
    constexpr u64 RADIX_MASK = static_cast<u64>(RADIX_SIZE - 1);
    std::array<std::uint32_t, RADIX_SIZE> count{};

    for (int pass = 0; pass < 4; ++pass) {
        const int shift = pass * RADIX_BITS;
        count.fill(0);

        for (const u64 v : values) {
            ++count[static_cast<std::size_t>((v >> shift) & RADIX_MASK)];
        }

        std::uint32_t acc = 0;
        for (int i = 0; i < RADIX_SIZE; ++i) {
            const std::uint32_t c = count[static_cast<std::size_t>(i)];
            count[static_cast<std::size_t>(i)] = acc;
            acc += c;
        }

        for (const u64 v : values) {
            tmp[count[static_cast<std::size_t>((v >> shift) & RADIX_MASK)]++] = v;
        }
        values.swap(tmp);
    }
}

u64 mod_pow_u64(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0) {
        if (e & 1ULL) {
            result = static_cast<u64>((__uint128_t)result * cur % mod);
        }
        cur = static_cast<u64>((__uint128_t)cur * cur % mod);
        e >>= 1ULL;
    }
    return result;
}

std::vector<int> unique_prime_factors(int x, const std::vector<int>& spf) {
    std::vector<int> fac;
    while (x > 1) {
        const int p = spf[static_cast<std::size_t>(x)];
        fac.push_back(p);
        while (x % p == 0) {
            x /= p;
        }
    }
    return fac;
}

int euler_phi(int x, const std::vector<int>& spf) {
    if (x == 1) {
        return 1;
    }
    int phi = x;
    int n = x;
    while (n > 1) {
        const int p = spf[static_cast<std::size_t>(n)];
        phi = phi / p * (p - 1);
        while (n % p == 0) {
            n /= p;
        }
    }
    return phi;
}

int multiplicative_order(int base, int mod, const std::vector<int>& spf) {
    if (mod == 1) {
        return 1;
    }
    int ord = euler_phi(mod, spf);
    std::vector<int> fac = unique_prime_factors(ord, spf);
    for (const int p : fac) {
        while (ord % p == 0) {
            const int cand = ord / p;
            if (mod_pow_u64(static_cast<u64>(base), static_cast<u64>(cand), static_cast<u64>(mod)) == 1ULL) {
                ord = cand;
            } else {
                break;
            }
        }
    }
    return ord;
}

int S_value(int n, const std::vector<int>& spf) {
    if (n == 1) {
        return 1;
    }

    int pre2 = 0;
    int tmp = n;
    while ((tmp & 1) == 0) {
        ++pre2;
        tmp >>= 1;
    }
    const int mod2 = tmp;
    const int per2 = multiplicative_order(2, mod2, spf);

    int pre3 = 0;
    tmp = n;
    while (tmp % 3 == 0) {
        ++pre3;
        tmp /= 3;
    }
    const int mod3 = tmp;
    const int per3 = multiplicative_order(3, mod3, spf);

    const int pre = std::max(pre2, pre3);
    const i64 period = std::lcm(static_cast<i64>(per2), static_cast<i64>(per3));
    const i64 needed = std::min<i64>(2LL * n + 1, pre + period);

    std::vector<u64> points;
    points.reserve(static_cast<std::size_t>(needed));

    int x = 1 % n;
    int y = 1 % n;
    for (i64 i = 0; i < needed; ++i) {
        const u64 key = (static_cast<u64>(static_cast<std::uint32_t>(x)) << 32) |
                        static_cast<u64>(static_cast<std::uint32_t>(y));
        points.push_back(key);
        x = static_cast<int>((2LL * x) % n);
        y = static_cast<int>((3LL * y) % n);
    }

    radix_sort_u64(points);
    points.erase(std::unique(points.begin(), points.end()), points.end());

    std::vector<int> tails;
    tails.reserve(256);
    for (const u64 key : points) {
        const int yy = static_cast<int>(key & 0xffffffffULL);
        auto it = std::upper_bound(tails.begin(), tails.end(), yy);
        if (it == tails.end()) {
            tails.push_back(yy);
        } else {
            *it = yy;
        }
    }

    return static_cast<int>(tails.size());
}

i64 solve_sum(int k_max, const std::vector<int>& spf) {
    i64 ans = 0;
    for (int k = 1; k <= k_max; ++k) {
        int n = 1;
        for (int t = 0; t < 5; ++t) {
            n *= k;
        }
        ans += S_value(n, spf);
    }
    return ans;
}

bool run_checkpoints(const std::vector<int>& spf) {
    if (S_value(22, spf) != 5) {
        std::cerr << "Checkpoint failed: S(22)\n";
        return false;
    }
    if (S_value(123, spf) != 14) {
        std::cerr << "Checkpoint failed: S(123)\n";
        return false;
    }
    if (S_value(10000, spf) != 48) {
        std::cerr << "Checkpoint failed: S(10000)\n";
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

    const int n_max = 24300000;  // 30^5
    const std::vector<int> spf = build_spf(n_max);

    if (options.run_checkpoints && !run_checkpoints(spf)) {
        return 2;
    }

    std::cout << solve_sum(options.k_max, spf) << '\n';
    return 0;
}
