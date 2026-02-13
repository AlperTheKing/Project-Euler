#include <cstdint>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using boost::multiprecision::cpp_int;

constexpr i64 MOD = 1'000'000'007LL;

i64 mod_pow(i64 base, i64 exp) {
    i64 result = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp & 1LL) {
            result = static_cast<i64>((__int128)result * base % MOD);
        }
        base = static_cast<i64>((__int128)base * base % MOD);
        exp >>= 1LL;
    }
    return result;
}

bool is_prime_small(int n) {
    if (n < 2) {
        return false;
    }
    for (int d = 2; static_cast<i64>(d) * d <= n; ++d) {
        if (n % d == 0) {
            return false;
        }
    }
    return true;
}

i64 sum_powers_mod(u64 n, int k) {
    if (n == 0) {
        return 0;
    }

    const int d = k + 2;
    std::vector<i64> y(static_cast<std::size_t>(d), 0);
    for (int i = 1; i < d; ++i) {
        y[static_cast<std::size_t>(i)] =
            (y[static_cast<std::size_t>(i - 1)] + mod_pow(i, k)) % MOD;
    }

    const i64 x = static_cast<i64>(n % static_cast<u64>(MOD));
    if (x < d) {
        return y[static_cast<std::size_t>(x)];
    }

    std::vector<i64> fact(static_cast<std::size_t>(d), 1);
    std::vector<i64> inv_fact(static_cast<std::size_t>(d), 1);
    for (int i = 1; i < d; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<i64>((__int128)fact[static_cast<std::size_t>(i - 1)] * i % MOD);
    }
    inv_fact[static_cast<std::size_t>(d - 1)] =
        mod_pow(fact[static_cast<std::size_t>(d - 1)], MOD - 2);
    for (int i = d - 1; i > 0; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<i64>((__int128)inv_fact[static_cast<std::size_t>(i)] * i % MOD);
    }

    std::vector<i64> prefix(static_cast<std::size_t>(d + 1), 1);
    std::vector<i64> suffix(static_cast<std::size_t>(d + 1), 1);

    for (int i = 0; i < d; ++i) {
        const i64 term = (x - i + MOD) % MOD;
        prefix[static_cast<std::size_t>(i + 1)] =
            static_cast<i64>((__int128)prefix[static_cast<std::size_t>(i)] * term % MOD);
    }
    for (int i = d - 1; i >= 0; --i) {
        const i64 term = (x - i + MOD) % MOD;
        suffix[static_cast<std::size_t>(i)] =
            static_cast<i64>((__int128)suffix[static_cast<std::size_t>(i + 1)] * term % MOD);
    }

    i64 answer = 0;
    for (int i = 0; i < d; ++i) {
        const i64 numerator =
            static_cast<i64>((__int128)prefix[static_cast<std::size_t>(i)] *
                                 suffix[static_cast<std::size_t>(i + 1)] %
                             MOD);

        i64 denominator = static_cast<i64>((__int128)inv_fact[static_cast<std::size_t>(i)] *
                                           inv_fact[static_cast<std::size_t>(d - 1 - i)] %
                                           MOD);
        if ((d - 1 - i) & 1) {
            denominator = (MOD - denominator) % MOD;
        }

        const i64 add = static_cast<i64>((__int128)y[static_cast<std::size_t>(i)] * numerator % MOD *
                                         denominator % MOD);
        answer += add;
        if (answer >= MOD) {
            answer -= MOD;
        }
    }

    return answer;
}

i64 solve_prefix(u64 n) {
    if (n == 0) {
        return 0;
    }

    i64 answer = 0;
    u128 primorial = 1;
    i64 primorial_mod = 1;
    i64 factorial_mod = 1;

    for (int m = 0;; ++m) {
        if (m >= 2 && is_prime_small(m)) {
            primorial *= static_cast<u128>(m);
            primorial_mod = static_cast<i64>((__int128)primorial_mod * m % MOD);
        }

        if (m > 0) {
            factorial_mod = static_cast<i64>((__int128)factorial_mod * m % MOD);
        }

        const u64 mu = static_cast<u64>(m);
        if (mu > n) {
            break;
        }
        const u64 remaining = n - mu;
        if (primorial > static_cast<u128>(remaining)) {
            break;
        }

        const u64 count = static_cast<u64>(static_cast<u128>(remaining) / primorial);

        i64 coefficient =
            static_cast<i64>((__int128)mod_pow(primorial_mod, m) *
                             mod_pow(factorial_mod, MOD - 2) %
                             MOD);
        if (m & 1) {
            coefficient = (MOD - coefficient) % MOD;
        }

        const i64 power_sum = sum_powers_mod(count, m);
        const i64 add = static_cast<i64>((__int128)coefficient * power_sum % MOD);

        answer += add;
        if (answer >= MOD) {
            answer -= MOD;
        }
    }

    return answer;
}

cpp_int direct_s(int n) {
    std::vector<cpp_int> factorial(static_cast<std::size_t>(n + 1), 1);
    for (int i = 1; i <= n; ++i) {
        factorial[static_cast<std::size_t>(i)] =
            factorial[static_cast<std::size_t>(i - 1)] * i;
    }

    cpp_int total = 0;
    for (int t = 1; t <= n; ++t) {
        const int m = n - t;
        cpp_int numerator = 1;
        for (int i = 0; i < m; ++i) {
            numerator *= t;
        }
        if (m & 1) {
            numerator = -numerator;
        }

        const cpp_int& denominator = factorial[static_cast<std::size_t>(m)];
        if (numerator % denominator == 0) {
            total += numerator / denominator;
        }
    }

    return total;
}

i64 cpp_int_mod(const cpp_int& value) {
    cpp_int residue = value % MOD;
    if (residue < 0) {
        residue += MOD;
    }
    return residue.convert_to<i64>();
}

bool run_checkpoints() {
    if (direct_s(1) != 1) {
        std::cerr << "Checkpoint failed: S(1) != 1" << '\n';
        return false;
    }

    if (direct_s(3) != -1) {
        std::cerr << "Checkpoint failed: S(3) != -1" << '\n';
        return false;
    }

    cpp_int prefix = 0;
    for (int n = 1; n <= 10; ++n) {
        prefix += direct_s(n);
    }
    if (prefix != 43) {
        std::cerr << "Checkpoint failed: sum_{n=1}^{10} S(n) != 43" << '\n';
        return false;
    }

    if (solve_prefix(10) != 43) {
        std::cerr << "Checkpoint failed: fast sum_{n=1}^{10} S(n) != 43" << '\n';
        return false;
    }

    for (int limit : {20, 30, 40}) {
        cpp_int direct_prefix = 0;
        for (int n = 1; n <= limit; ++n) {
            direct_prefix += direct_s(n);
        }

        const i64 fast = solve_prefix(static_cast<u64>(limit));
        const i64 slow = cpp_int_mod(direct_prefix);
        if (fast != slow) {
            std::cerr << "Checkpoint failed for N=" << limit << ": fast=" << fast
                      << ", slow=" << slow << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 target = 1'000'000'000'000'000'000ULL;
    std::cout << solve_prefix(target) << '\n';
    return 0;
}
