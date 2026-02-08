#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int limit = 2'000'000;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
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

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

u64 tonelli_shanks(const u64 n, const u64 p) {
    if (n == 0ULL) {
        return 0ULL;
    }
    if (p == 2ULL) {
        return n;
    }
    if (pow_mod(n, (p - 1ULL) / 2ULL, p) != 1ULL) {
        return 0ULL;
    }
    if ((p & 3ULL) == 3ULL) {
        return pow_mod(n, (p + 1ULL) / 4ULL, p);
    }

    u64 q = p - 1ULL;
    int s = 0;
    while ((q & 1ULL) == 0ULL) {
        q >>= 1ULL;
        ++s;
    }

    u64 z = 2ULL;
    while (pow_mod(z, (p - 1ULL) / 2ULL, p) != p - 1ULL) {
        ++z;
    }

    u64 c = pow_mod(z, q, p);
    u64 x = pow_mod(n, (q + 1ULL) / 2ULL, p);
    u64 t = pow_mod(n, q, p);
    int m = s;

    while (t != 1ULL) {
        int i = 1;
        u64 t2i = mul_mod(t, t, p);
        while (i < m && t2i != 1ULL) {
            t2i = mul_mod(t2i, t2i, p);
            ++i;
        }
        const u64 b = pow_mod(c, 1ULL << (m - i - 1), p);
        x = mul_mod(x, b, p);
        t = mul_mod(t, mul_mod(b, b, p), p);
        c = mul_mod(b, b, p);
        m = i;
    }
    return x;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int d = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + d));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 solve(const int limit) {
    std::vector<int> lpf_small(static_cast<std::size_t>(limit + 2), 0);
    std::vector<int> primes;
    for (int i = 2; i <= limit + 1; ++i) {
        if (lpf_small[static_cast<std::size_t>(i)] == 0) {
            for (int j = i; j <= limit + 1; j += i) {
                lpf_small[static_cast<std::size_t>(j)] = i;
            }
            if (i <= limit) {
                primes.push_back(i);
            }
        }
    }

    std::vector<u64> rem(static_cast<std::size_t>(limit + 1), 0ULL);
    std::vector<u64> lpf_q(static_cast<std::size_t>(limit + 1), 1ULL);
    for (int k = 1; k <= limit; ++k) {
        const u64 kk = static_cast<u64>(k);
        rem[static_cast<std::size_t>(k)] = kk * kk - kk + 1ULL;
    }

    for (int p : primes) {
        if (p == 2) {
            continue;
        }
        std::vector<int> roots;
        if (p == 3) {
            roots.push_back(2);
        } else {
            const u64 mod = static_cast<u64>(p);
            const u64 neg_three = (mod + mod - 3ULL) % mod;
            if (pow_mod(neg_three, (mod - 1ULL) / 2ULL, mod) != 1ULL) {
                continue;
            }
            const u64 sq = tonelli_shanks(neg_three, mod);
            const u64 inv2 = (mod + 1ULL) / 2ULL;
            const int r1 = static_cast<int>(mul_mod((1ULL + sq) % mod, inv2, mod));
            const int r2 = static_cast<int>(mul_mod((1ULL + mod - sq) % mod, inv2, mod));
            roots.push_back(r1);
            if (r2 != r1) {
                roots.push_back(r2);
            }
        }

        for (int root : roots) {
            int k = root;
            if (k == 0) {
                k += p;
            }
            for (; k <= limit; k += p) {
                u64& value = rem[static_cast<std::size_t>(k)];
                if (value % static_cast<u64>(p) != 0ULL) {
                    continue;
                }
                while (value % static_cast<u64>(p) == 0ULL) {
                    value /= static_cast<u64>(p);
                }
                lpf_q[static_cast<std::size_t>(k)] = static_cast<u64>(p);
            }
        }
    }

    u128 total = 0;
    for (int k = 1; k <= limit; ++k) {
        const u64 rem_val = rem[static_cast<std::size_t>(k)];
        if (rem_val > 1ULL) {
            lpf_q[static_cast<std::size_t>(k)] = std::max(lpf_q[static_cast<std::size_t>(k)], rem_val);
        }
        const u64 lpf_kp1 = static_cast<u64>(lpf_small[static_cast<std::size_t>(k + 1)]);
        const u64 largest = std::max(lpf_kp1, lpf_q[static_cast<std::size_t>(k)]);
        total += static_cast<u128>(largest - 1ULL);
    }
    return total;
}

u64 brute_f_of_k(const int k) {
    u64 x = 1ULL;
    u64 y = static_cast<u64>(k);
    while (y != 1ULL) {
        ++x;
        --y;
        const u64 g = std::gcd(x, y);
        x /= g;
        y /= g;
    }
    return x;
}

u128 brute_sum_k_cubed(const int limit) {
    u128 s = 0;
    for (int k = 1; k <= limit; ++k) {
        const int kc = k * k * k;
        s += static_cast<u128>(brute_f_of_k(kc));
    }
    return s;
}

bool run_checkpoints() {
    if (brute_f_of_k(20) != 6ULL) {
        std::cerr << "Checkpoint failed: f(20)=6" << '\n';
        return false;
    }
    if (to_string_u128(solve(100)) != "118937") {
        std::cerr << "Checkpoint failed: sum_{k<=100} f(k^3)=118937" << '\n';
        return false;
    }
    if (solve(300) != brute_sum_k_cubed(300)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for limit=300" << '\n';
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
    std::cout << to_string_u128(solve(options.limit)) << '\n';
    return 0;
}
