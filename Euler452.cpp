#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u32 kMod = 1'234'567'891U;

struct Options {
    u32 m = 1'000'000'000U;
    u32 n = 1'000'000'000U;
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
        if (parse_u32_after_prefix(arg, "--m=", options.m)) {
            continue;
        }
        if (parse_u32_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.m >= 1U;
}

u32 choose_small_mod(const u32 n, const u32 r) {
    if (r > n) {
        return 0U;
    }
    if (r == 0U) {
        return 1U;
    }

    std::vector<u64> numerators(r);
    for (u32 i = 0U; i < r; ++i) {
        numerators[static_cast<std::size_t>(i)] = static_cast<u64>(n - r + 1U + i);
    }

    for (u32 d = 2U; d <= r; ++d) {
        u64 x = d;
        for (u32 i = 0U; i < r && x > 1ULL; ++i) {
            u64& v = numerators[static_cast<std::size_t>(i)];
            const u64 g = std::gcd(v, x);
            if (g > 1ULL) {
                v /= g;
                x /= g;
            }
        }
    }

    u64 result = 1ULL;
    for (const u64 v : numerators) {
        result = (result * (v % kMod)) % kMod;
    }
    return static_cast<u32>(result);
}

u32 solve(const u32 m, const u32 n_target) {
    int degree = 0;
    while ((1ULL << (degree + 1)) <= m) {
        ++degree;
    }

    std::vector<u32> values;
    values.reserve(static_cast<std::size_t>(2.0L * std::sqrt(static_cast<long double>(m)) + 8.0L));
    for (u32 i = 1U; i <= m;) {
        const u32 q = m / i;
        const u32 j = m / q;
        values.push_back(q);
        i = j + 1U;
    }

    const u32 state_count = static_cast<u32>(values.size());

    u32 sqrt_m = static_cast<u32>(std::sqrt(static_cast<long double>(m)));
    while ((static_cast<u64>(sqrt_m) + 1ULL) * (static_cast<u64>(sqrt_m) + 1ULL) <= m) {
        ++sqrt_m;
    }
    while (static_cast<u64>(sqrt_m) * static_cast<u64>(sqrt_m) > m) {
        --sqrt_m;
    }

    const auto index_of = [&](const u32 q) -> u32 {
        if (q <= sqrt_m) {
            return state_count - q;
        }
        return m / q - 1U;
    };

    std::vector<u32> offset(static_cast<std::size_t>(state_count) + 1U, 0U);
    std::vector<u32> q_index;
    std::vector<u32> count_mod;

    for (u32 idx = 0U; idx < state_count; ++idx) {
        const u32 v = values[static_cast<std::size_t>(idx)];
        offset[static_cast<std::size_t>(idx)] = static_cast<u32>(q_index.size());

        for (u32 i = 1U; i <= v;) {
            const u32 q = v / i;
            const u32 j = v / q;
            q_index.push_back(index_of(q));
            count_mod.push_back((j - i + 1U) % kMod);
            i = j + 1U;
        }
    }
    offset[static_cast<std::size_t>(state_count)] = static_cast<u32>(q_index.size());

    std::vector<u32> prev(state_count, 1U);   // F(v,0)=1
    std::vector<u32> curr(state_count, 0U);
    std::vector<u32> poly_values(static_cast<std::size_t>(degree) + 1U, 0U);
    poly_values[0] = 1U;

    for (int t = 1; t <= degree; ++t) {
        for (u32 idx = 0U; idx < state_count; ++idx) {
            const u32 begin = offset[static_cast<std::size_t>(idx)];
            const u32 end = offset[static_cast<std::size_t>(idx) + 1U];

            u128 sum = 0;
            int chunk = 0;
            for (u32 p = begin; p < end; ++p) {
                sum += static_cast<u128>(count_mod[static_cast<std::size_t>(p)]) *
                       prev[static_cast<std::size_t>(q_index[static_cast<std::size_t>(p)])];
                ++chunk;
                if (chunk == 64) {
                    sum %= kMod;
                    chunk = 0;
                }
            }
            curr[static_cast<std::size_t>(idx)] = static_cast<u32>(sum % kMod);
        }

        prev.swap(curr);
        poly_values[static_cast<std::size_t>(t)] = prev[0];
    }

    std::vector<u32> diff = poly_values;
    std::vector<u32> coeff(static_cast<std::size_t>(degree) + 1U, 0U);
    for (int r = 0; r <= degree; ++r) {
        coeff[static_cast<std::size_t>(r)] = diff[0];
        for (int i = 0; i < degree - r; ++i) {
            u32 v = diff[static_cast<std::size_t>(i + 1)] + kMod - diff[static_cast<std::size_t>(i)];
            if (v >= kMod) {
                v -= kMod;
            }
            diff[static_cast<std::size_t>(i)] = v;
        }
    }

    u64 answer = 0ULL;
    for (int r = 0; r <= degree; ++r) {
        const u32 c = coeff[static_cast<std::size_t>(r)];
        const u32 binom = choose_small_mod(n_target, static_cast<u32>(r));
        answer += static_cast<u64>(c) * binom;
        answer %= kMod;
    }

    return static_cast<u32>(answer);
}

bool run_checkpoints() {
    if (solve(10U, 10U) != 571U) {
        std::cerr << "Checkpoint failed: F(10,10)\n";
        return false;
    }
    if (solve(1'000'000U, 1'000'000U) != 252'903'833U) {
        std::cerr << "Checkpoint failed: F(10^6,10^6)\n";
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

    std::cout << solve(options.m, options.n) << '\n';
    return 0;
}
