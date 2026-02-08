#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
    u64 a = 1'801'088'541ULL;                    // 21^7
    u64 b = 558'545'864'083'284'007ULL;          // 7^21
    u64 c = 35'831'808ULL;                       // 12^7
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
    u64 parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
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
        if (parse_u64_after_prefix(arg, "--a=", options.a) ||
            parse_u64_after_prefix(arg, "--b=", options.b) ||
            parse_u64_after_prefix(arg, "--c=", options.c)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.a > options.c && options.b >= options.a;
}

u64 F_closed(const u64 n, const u64 a, const u64 b, const u64 c) {
    if (n > b) {
        return n - c;
    }
    const u64 d = a - c;
    const u64 q = (b - n) / a;
    return n + 4ULL * d + (a + 3ULL * d) * q;
}

u64 S_closed_mod_1e9(const u64 a, const u64 b, const u64 c) {
    constexpr u64 MOD = 1'000'000'000ULL;
    const u64 d = a - c;
    const u64 q = b / a;
    const u64 r = b % a;

    const auto mod_mul = [](u64 x, u64 y) -> u64 {
        return static_cast<u64>((static_cast<__uint128_t>(x) * y) % 1'000'000'000ULL);
    };

    u64 sum_n = static_cast<u64>((static_cast<__uint128_t>(b % MOD) * ((b + 1ULL) % MOD) / 2ULL) % MOD);
    if ((b & 1ULL) == 0ULL) {
        sum_n = mod_mul((b / 2ULL) % MOD, (b + 1ULL) % MOD);
    } else {
        sum_n = mod_mul(b % MOD, ((b + 1ULL) / 2ULL) % MOD);
    }

    u64 term_linear = mod_mul((b + 1ULL) % MOD, (4ULL * (d % MOD)) % MOD);

    i128 sum_floor =
        static_cast<i128>(a) * static_cast<i128>(q) * static_cast<i128>(q - 1ULL) / 2 +
        static_cast<i128>(r + 1ULL) * static_cast<i128>(q);
    const u64 sum_floor_mod = static_cast<u64>(sum_floor % MOD);
    u64 term_floor = mod_mul((a + 3ULL * d) % MOD, sum_floor_mod);

    return (sum_n + term_linear + term_floor) % MOD;
}

u64 S_closed_exact_small(const u64 a, const u64 b, const u64 c) {
    const u64 d = a - c;
    const u64 q = b / a;
    const u64 r = b % a;
    const i128 sum_floor =
        static_cast<i128>(a) * static_cast<i128>(q) * static_cast<i128>(q - 1ULL) / 2 +
        static_cast<i128>(r + 1ULL) * static_cast<i128>(q);
    const i128 result =
        static_cast<i128>(b) * static_cast<i128>(b + 1ULL) / 2 +
        static_cast<i128>(b + 1ULL) * static_cast<i128>(4ULL * d) +
        static_cast<i128>(a + 3ULL * d) * sum_floor;
    return static_cast<u64>(result);
}

u64 F_bruteforce(const u64 n, const u64 a, const u64 b, const u64 c) {
    std::unordered_map<u64, u64> memo;
    std::function<u64(u64)> rec = [&](u64 x) -> u64 {
        auto it = memo.find(x);
        if (it != memo.end()) {
            return it->second;
        }
        u64 out = 0ULL;
        if (x > b) {
            out = x - c;
        } else {
            out = rec(a + rec(a + rec(a + rec(a + x))));
        }
        memo.emplace(x, out);
        return out;
    };
    return rec(n);
}

u64 S_bruteforce(const u64 a, const u64 b, const u64 c) {
    u64 sum = 0ULL;
    for (u64 n = 0ULL; n <= b; ++n) {
        sum += F_bruteforce(n, a, b, c);
    }
    return sum;
}

bool run_checkpoints() {
    if (F_closed(0, 50, 2000, 40) != 3240ULL) {
        std::cerr << "Checkpoint failed: F(0) for sample parameters" << '\n';
        return false;
    }
    if (F_closed(2000, 50, 2000, 40) != 2040ULL) {
        std::cerr << "Checkpoint failed: F(2000) for sample parameters" << '\n';
        return false;
    }
    if (S_closed_exact_small(50, 2000, 40) != 5'204'240ULL) {
        std::cerr << "Checkpoint failed: S(50,2000,40)" << '\n';
        return false;
    }
    for (u64 a = 4; a <= 9; ++a) {
        for (u64 c = 1; c < a; ++c) {
            const u64 b = 60 + a;
            for (u64 n = 0; n <= b; ++n) {
                if (F_closed(n, a, b, c) != F_bruteforce(n, a, b, c)) {
                    std::cerr << "Checkpoint failed: brute-force cross-check for F(n)" << '\n';
                    return false;
                }
            }
            if (S_closed_exact_small(a, b, c) != S_bruteforce(a, b, c)) {
                std::cerr << "Checkpoint failed: brute-force cross-check for S" << '\n';
                return false;
            }
        }
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
    std::cout << S_closed_mod_1e9(options.a, options.b, options.c) << '\n';
    return 0;
}
