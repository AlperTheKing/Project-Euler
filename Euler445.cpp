#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kMod = 1'000'000'007U;

struct Options {
    u32 n = 10'000'000U;
    bool run_checkpoints = true;
};

bool parse_u32_after_prefix(const std::string& arg, const std::string& prefix, u32& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        const unsigned long long parsed = std::stoull(tail);
        if (parsed > static_cast<unsigned long long>(std::numeric_limits<u32>::max())) {
            return false;
        }
        out = static_cast<u32>(parsed);
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
        if (parse_u32_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 2U;
}

u32 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    base %= kMod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = (result * base) % kMod;
        }
        base = (base * base) % kMod;
        exp >>= 1ULL;
    }
    return static_cast<u32>(result);
}

class TermProductTracker {
public:
    void add(u32 value) {
        if (value == 0U) {
            ++zero_terms_;
            return;
        }
        product_non_zero_ = (product_non_zero_ * value) % kMod;
    }

    void remove(u32 value) {
        if (value == 0U) {
            --zero_terms_;
            return;
        }
        product_non_zero_ = (product_non_zero_ * inverse(value)) % kMod;
    }

    u32 value() const {
        return (zero_terms_ > 0) ? 0U : static_cast<u32>(product_non_zero_);
    }

private:
    u64 product_non_zero_ = 1ULL;
    int zero_terms_ = 0;
    std::unordered_map<u32, u32> inverse_cache_;

    u32 inverse(const u32 x) {
        const auto it = inverse_cache_.find(x);
        if (it != inverse_cache_.end()) {
            return it->second;
        }
        const u32 inv = mod_pow(x, static_cast<u64>(kMod) - 2ULL);
        inverse_cache_.emplace(x, inv);
        return inv;
    }
};

u32 solve(const u32 n) {
    std::vector<u32> spf(static_cast<std::size_t>(n) + 1U, 0U);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(n / 10U));

    for (u32 i = 2U; i <= n; ++i) {
        if (spf[i] == 0U) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (const u32 p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > n || p > spf[i]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }

    std::vector<u32> inv(static_cast<std::size_t>(n) + 1U, 0U);
    inv[1] = 1U;
    for (u32 i = 2U; i <= n; ++i) {
        inv[i] = static_cast<u32>(kMod - (static_cast<u64>(kMod / i) * inv[kMod % i]) % kMod);
    }

    std::vector<int> exponent(static_cast<std::size_t>(n) + 1U, 0);
    std::vector<u32> power_mod(static_cast<std::size_t>(n) + 1U, 1U);
    TermProductTracker tracker;

    const auto adjust_prime = [&](const u32 p, const int delta) {
        int old_exp = exponent[p];
        u32 p_pow = power_mod[p];

        if (old_exp > 0) {
            tracker.remove((p_pow + 1U) % kMod);
        }

        if (delta > 0) {
            for (int i = 0; i < delta; ++i) {
                p_pow = static_cast<u32>((static_cast<u64>(p_pow) * p) % kMod);
            }
        } else if (delta < 0) {
            const u32 inv_p = inv[p];
            for (int i = 0; i < -delta; ++i) {
                p_pow = static_cast<u32>((static_cast<u64>(p_pow) * inv_p) % kMod);
            }
        }

        old_exp += delta;
        exponent[p] = old_exp;
        power_mod[p] = p_pow;

        if (old_exp > 0) {
            tracker.add((p_pow + 1U) % kMod);
        }
    };

    const auto apply_factorization = [&](u32 x, const int sign) {
        while (x > 1U) {
            const u32 p = spf[x];
            int cnt = 0;
            while (x % p == 0U) {
                x /= p;
                ++cnt;
            }
            adjust_prime(p, sign * cnt);
        }
    };

    u64 binom_mod = 1ULL;  // C(n,0)
    u64 answer = 0ULL;
    const u32 half = n / 2U;

    for (u32 k = 1U; k <= half; ++k) {
        apply_factorization(n - k + 1U, +1);
        apply_factorization(k, -1);

        binom_mod = (binom_mod * (n - k + 1ULL)) % kMod;
        binom_mod = (binom_mod * inv[k]) % kMod;

        const u32 q_mod = tracker.value();
        const u32 r_mod = static_cast<u32>((q_mod + kMod - binom_mod) % kMod);
        const u32 weight = (k == n - k) ? 1U : 2U;
        answer += static_cast<u64>(weight) * r_mod;
        answer %= kMod;
    }

    return static_cast<u32>(answer);
}

bool run_checkpoints() {
    if (solve(10U) != 1118U) {
        std::cerr << "Checkpoint failed: N=10\n";
        return false;
    }
    if (solve(100'000U) != 628'701'600U) {
        std::cerr << "Checkpoint failed: N=100000\n";
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
