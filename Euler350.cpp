#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 MOD = 104060401;  // 101^4

i64 mod_add(const i64 a, const i64 b) {
    i64 x = a + b;
    if (x >= MOD) {
        x -= MOD;
    }
    return x;
}

i64 mod_sub(const i64 a, const i64 b) {
    i64 x = a - b;
    if (x < 0) {
        x += MOD;
    }
    return x;
}

i64 mod_mul(const i64 a, const i64 b) {
    return static_cast<i64>((static_cast<__int128>(a) * static_cast<__int128>(b)) % MOD);
}

i64 mod_pow(i64 base, u64 exp) {
    base %= MOD;
    i64 result = 1;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

std::vector<int> build_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1));
    for (int i = 0; i <= n; ++i) {
        spf[static_cast<std::size_t>(i)] = i;
    }

    for (int i = 2; static_cast<i64>(i) * i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            if (spf[static_cast<std::size_t>(j)] == j) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }

    return spf;
}

i64 solve(const u64 G, const u64 L, const u64 N) {
    const int K = static_cast<int>(L / G);
    const std::vector<int> spf = build_spf(K);

    int max_exp = 0;
    {
        int x = K;
        while (x > 0) {
            x /= 2;
            ++max_exp;
        }
    }

    std::vector<i64> pows(static_cast<std::size_t>(max_exp + 2), 0);
    for (int b = 0; b <= max_exp + 1; ++b) {
        if (b == 0) {
            pows[static_cast<std::size_t>(b)] = 0;  // N > 0 here
        } else {
            pows[static_cast<std::size_t>(b)] = mod_pow(b, N);
        }
    }

    std::vector<i64> local(static_cast<std::size_t>(max_exp + 1), 0);
    for (int a = 1; a <= max_exp; ++a) {
        // (a+1)^N - 2*a^N + (a-1)^N
        i64 value = mod_sub(pows[static_cast<std::size_t>(a + 1)], mod_mul(2, pows[static_cast<std::size_t>(a)]));
        value = mod_add(value, pows[static_cast<std::size_t>(a - 1)]);
        local[static_cast<std::size_t>(a)] = value;
    }

    std::vector<i64> h(static_cast<std::size_t>(K + 1), 0);
    h[1] = 1;

    for (int n = 2; n <= K; ++n) {
        int x = n;
        i64 value = 1;

        while (x > 1) {
            const int p = spf[static_cast<std::size_t>(x)];
            int e = 0;
            while (x % p == 0) {
                x /= p;
                ++e;
            }
            value = mod_mul(value, local[static_cast<std::size_t>(e)]);
        }

        h[static_cast<std::size_t>(n)] = value;
    }

    i64 answer = 0;
    for (int m = 1; m <= K; ++m) {
        const u64 weight = L / static_cast<u64>(m) - G + 1ULL;
        const i64 w = static_cast<i64>(weight % static_cast<u64>(MOD));
        answer = mod_add(answer, mod_mul(h[static_cast<std::size_t>(m)], w));
    }

    return answer;
}

bool run_checkpoints() {
    if (solve(10ULL, 100ULL, 1ULL) != 91) {
        std::cerr << "Checkpoint failed: f(10,100,1)\n";
        return false;
    }
    if (solve(10ULL, 100ULL, 2ULL) != 327) {
        std::cerr << "Checkpoint failed: f(10,100,2)\n";
        return false;
    }
    if (solve(10ULL, 100ULL, 3ULL) != 1135) {
        std::cerr << "Checkpoint failed: f(10,100,3)\n";
        return false;
    }
    if (solve(10ULL, 100ULL, 1000ULL) != 3286053) {
        std::cerr << "Checkpoint failed: f(10,100,1000) mod 101^4\n";
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

    const i64 answer = solve(1000000ULL, 1000000000000ULL, 1000000000000000000ULL);
    std::cout << answer << '\n';
    return 0;
}
