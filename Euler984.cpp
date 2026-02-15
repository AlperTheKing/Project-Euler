#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

static constexpr i64 MOD = 1000000007LL;
static constexpr int K = 14;

static constexpr std::array<i64, K> INIT = {
    1LL,
    4LL,
    9LL,
    92LL,
    903LL,
    4411LL,
    14959LL,
    41083LL,
    98200LL,
    212418LL,
    425756LL,
    803074LL,
    1441065LL,
    2479669LL
};

static constexpr std::array<i64, K> REC = {
    7LL,
    MOD - 19LL,
    21LL,
    6LL,
    MOD - 42LL,
    42LL,
    MOD - 6LL,
    MOD - 21LL,
    19LL,
    MOD - 7LL,
    1LL,
    0LL,
    0LL,
    0LL
};

static constexpr std::array<i64, 11> REC_SIGNED = {
    7LL,
    -19LL,
    21LL,
    6LL,
    -42LL,
    42LL,
    -6LL,
    -21LL,
    19LL,
    -7LL,
    1LL
};

static i64 mod_pow(i64 base, u64 exp) {
    i64 result = 1;
    base %= MOD;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = static_cast<i128>(result) * base % MOD;
        }
        base = static_cast<i128>(base) * base % MOD;
        exp >>= 1ULL;
    }
    return result;
}

static i64 mod_inv(i64 x) {
    return mod_pow((x % MOD + MOD) % MOD, MOD - 2);
}

static i64 f_even_closed_form_mod(const u64 n) {
    const i64 x = static_cast<i64>(n % static_cast<u64>(MOD));
    std::array<i64, 9> p{};
    p[0] = 1;
    for (int i = 1; i <= 8; ++i) {
        p[i] = static_cast<i128>(p[i - 1]) * x % MOD;
    }

    const i64 inv40320 = mod_inv(40320);
    const i64 inv3360 = mod_inv(3360);
    const i64 inv1440 = mod_inv(1440);
    const i64 inv320 = mod_inv(320);
    const i64 inv240 = mod_inv(240);
    const i64 inv420 = mod_inv(420);
    const i64 inv140 = mod_inv(140);

    i64 ans = 0;
    ans = (ans + static_cast<i128>(31) * p[8] % MOD * inv40320) % MOD;
    ans = (ans + static_cast<i128>(31) * p[7] % MOD * inv3360) % MOD;
    ans = (ans + static_cast<i128>(67) * p[6] % MOD * inv1440) % MOD;
    ans = (ans + static_cast<i128>(41) * p[5] % MOD * inv320) % MOD;
    ans = (ans + static_cast<i128>(313) * p[4] % MOD * inv1440) % MOD;
    ans = (ans + static_cast<i128>(MOD - 5699) * p[3] % MOD * inv240) % MOD;
    ans = (ans + static_cast<i128>(16049) * p[2] % MOD * inv420) % MOD;
    ans = (ans + static_cast<i128>(29413) * p[1] % MOD * inv140) % MOD;
    ans = (ans + MOD - 419) % MOD;
    return ans;
}

static std::array<i64, K> combine(const std::array<i64, K>& a,
                                   const std::array<i64, K>& b) {
    std::array<i64, 2 * K> prod{};

    for (int i = 0; i < K; ++i) {
        if (a[i] == 0) continue;
        for (int j = 0; j < K; ++j) {
            if (b[j] == 0) continue;
            prod[i + j] = (prod[i + j] + static_cast<i128>(a[i]) * b[j]) % MOD;
        }
    }

    for (int i = 2 * K - 2; i >= K; --i) {
        if (prod[i] == 0) continue;
        for (int j = 1; j <= K; ++j) {
            prod[i - j] = (prod[i - j] + static_cast<i128>(prod[i]) * REC[j - 1]) % MOD;
        }
    }

    std::array<i64, K> out{};
    for (int i = 0; i < K; ++i) out[i] = prod[i];
    return out;
}

static i64 f_mod(const u64 n) {
    if (n == 0ULL) {
        return 0;
    }
    if (n > 3ULL && (n & 1ULL) == 0ULL) {
        return f_even_closed_form_mod(n);
    }
    if (n <= static_cast<u64>(K)) {
        return INIT[static_cast<std::size_t>(n - 1)] % MOD;
    }

    u64 exp = n - 1ULL;
    std::array<i64, K> poly{};
    std::array<i64, K> step{};
    poly[0] = 1;
    step[1] = 1;

    while (exp > 0ULL) {
        if (exp & 1ULL) {
            poly = combine(poly, step);
        }
        step = combine(step, step);
        exp >>= 1ULL;
    }

    i64 ans = 0;
    for (int i = 0; i < K; ++i) {
        ans = (ans + static_cast<i128>(poly[i]) * INIT[i]) % MOD;
    }
    return ans;
}

static i64 f_100_exact() {
    std::vector<i128> f(101, 0);
    for (int n = 1; n <= K; ++n) {
        f[n] = INIT[static_cast<std::size_t>(n - 1)];
    }

    // This recurrence is valid from n >= 15.
    for (int n = 15; n <= 100; ++n) {
        i128 cur = 0;
        for (int j = 1; j <= 11; ++j) {
            cur += static_cast<i128>(REC_SIGNED[static_cast<std::size_t>(j - 1)]) * f[n - j];
        }
        f[n] = cur;
    }

    return static_cast<i64>(f[100]);
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (f_mod(3ULL) != 9LL ||
        f_mod(5ULL) != 903LL ||
        f_100_exact() != 8658918531876LL ||
        f_mod(10000ULL) != 377956308LL) {
        std::cerr << "Validation failed\n";
        return 1;
    }

    const u64 target = 1000000000000000000ULL;
    std::cout << f_mod(target) << '\n';
    return 0;
}
