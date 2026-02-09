#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 100'000'000ULL;

struct Options {
    u64 n = 1'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n == 0ULL) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

struct Precomp {
    // block k stores f(2^k + 1) ... f(2^(k+1)), length 2^k
    std::array<u128, 64> block_sum{};
    std::array<u128, 64> block_pref_3q{};
    std::array<u128, 64> pow_prefix{};
    std::array<std::vector<u64>, 5> base_blocks{};
    std::array<std::vector<u128>, 5> base_pref{};
};

u128 prefix_R(const u64 r, const u64 m) {
    if (m == 0ULL) {
        return 0;
    }
    if (m == 1ULL) {
        return static_cast<u128>(4ULL * r + 2ULL);
    }
    const u64 t = m - 1ULL;
    return static_cast<u128>(4ULL * r + 2ULL) + static_cast<u128>(t) *
                                                      static_cast<u128>(2ULL * r - m + 2ULL);
}

u128 sum_R(const u64 r) { return static_cast<u128>(r) * static_cast<u128>(r + 5ULL); }

u128 prefix_S(const u64 m) { return static_cast<u128>(m) * static_cast<u128>(m + 1ULL); }

u128 sum_S(const u64 r) {
    const u128 rr = static_cast<u128>(r) * static_cast<u128>(r);
    return 4ULL * rr + 2ULL * static_cast<u128>(r);
}

u128 prefix_U(const u64 r, const u64 m) {
    if (m == 0ULL) {
        return 0;
    }
    if (m == 1ULL) {
        return static_cast<u128>(6ULL * r + 3ULL);
    }
    const u64 t = m - 1ULL;
    return static_cast<u128>(6ULL * r + 3ULL) +
           (static_cast<u128>(t) * static_cast<u128>(4ULL * r - m + 6ULL)) / 2ULL;
}

u128 sum_U(const u64 r) {
    const u128 rr = static_cast<u128>(r) * static_cast<u128>(r);
    return 2ULL * rr + 11ULL * static_cast<u128>(r);
}

Precomp build_precomp() {
    Precomp pc;
    pc.base_blocks[0] = {2ULL};
    pc.base_blocks[1] = {2ULL, 4ULL};
    pc.base_blocks[2] = {3ULL, 6ULL, 2ULL, 6ULL};
    pc.base_blocks[3] = {3ULL, 8ULL, 2ULL, 6ULL, 2ULL, 4ULL, 9ULL, 4ULL};
    pc.base_blocks[4] = {3ULL, 8ULL, 2ULL, 6ULL, 2ULL, 4ULL, 10ULL, 4ULL,
                         2ULL, 4ULL, 6ULL, 8ULL, 15ULL, 6ULL, 5ULL, 4ULL};

    for (int k = 0; k <= 4; ++k) {
        const auto& blk = pc.base_blocks[static_cast<std::size_t>(k)];
        auto& pref = pc.base_pref[static_cast<std::size_t>(k)];
        pref.assign(blk.size() + 1U, 0);
        for (std::size_t i = 0; i < blk.size(); ++i) {
            pref[i + 1U] = pref[i] + static_cast<u128>(blk[i]);
        }
        pc.block_sum[static_cast<std::size_t>(k)] = pref.back();

        const std::size_t first_three_quarters = (3ULL * (1ULL << k)) / 4ULL;
        pc.block_pref_3q[static_cast<std::size_t>(k)] = pref[first_three_quarters];
    }

    for (int k = 5; k < 64; ++k) {
        const u64 r = (1ULL << (k - 3));
        const u128 rr = static_cast<u128>(r) * static_cast<u128>(r);
        pc.block_pref_3q[static_cast<std::size_t>(k)] =
            pc.block_pref_3q[static_cast<std::size_t>(k - 1)] + 5ULL * rr + 7ULL * r;
        pc.block_sum[static_cast<std::size_t>(k)] =
            pc.block_pref_3q[static_cast<std::size_t>(k)] + 2ULL * rr + 11ULL * r;
    }

    pc.pow_prefix[0] = 1;  // sum_{n<=1} f(n)
    for (int k = 1; k < 64; ++k) {
        pc.pow_prefix[static_cast<std::size_t>(k)] =
            pc.pow_prefix[static_cast<std::size_t>(k - 1)] +
            pc.block_sum[static_cast<std::size_t>(k - 1)];
    }
    return pc;
}

u128 block_prefix_sum(const Precomp& pc, const int k, const u64 m) {
    if (m == 0ULL) {
        return 0;
    }
    const u64 len = (1ULL << k);
    if (m >= len) {
        return pc.block_sum[static_cast<std::size_t>(k)];
    }

    if (k <= 4) {
        return pc.base_pref[static_cast<std::size_t>(k)][static_cast<std::size_t>(m)];
    }

    const u64 q = 3ULL * (1ULL << (k - 3));
    const u64 r = (1ULL << (k - 3));
    const u64 s = (1ULL << (k - 2));

    if (m <= q) {
        return block_prefix_sum(pc, k - 1, m);
    }
    if (m <= q + r) {
        return pc.block_pref_3q[static_cast<std::size_t>(k - 1)] + prefix_R(r, m - q);
    }
    if (m <= q + r + s) {
        return pc.block_pref_3q[static_cast<std::size_t>(k - 1)] + sum_R(r) +
               prefix_S(m - q - r);
    }
    return pc.block_pref_3q[static_cast<std::size_t>(k - 1)] + sum_R(r) + sum_S(r) +
           prefix_U(r, m - q - r - s);
}

u128 prefix_sum(const Precomp& pc, const u64 n) {
    if (n == 0ULL) {
        return 0;
    }
    if (n == 1ULL) {
        return 1;
    }
    const int k = 63 - __builtin_clzll(n);
    const u64 p = (1ULL << k);
    u128 out = pc.pow_prefix[static_cast<std::size_t>(k)];
    if (n > p) {
        out += block_prefix_sum(pc, k, n - p);
    }
    return out;
}

u64 f_value(const Precomp& pc, const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    return static_cast<u64>(prefix_sum(pc, n) - prefix_sum(pc, n - 1ULL));
}

bool run_checkpoints(const Precomp& pc) {
    if (f_value(pc, 1ULL) != 1ULL) {
        std::cerr << "Checkpoint failed: f(1)\n";
        return false;
    }
    if (f_value(pc, 15ULL) != 9ULL) {
        std::cerr << "Checkpoint failed: f(15)\n";
        return false;
    }
    if (f_value(pc, 20ULL) != 6ULL) {
        std::cerr << "Checkpoint failed: f(20)\n";
        return false;
    }
    if (f_value(pc, 500ULL) != 16ULL) {
        std::cerr << "Checkpoint failed: f(500)\n";
        return false;
    }
    if (prefix_sum(pc, 20ULL) != static_cast<u128>(83ULL)) {
        std::cerr << "Checkpoint failed: sum f(N), 1<=N<=20\n";
        return false;
    }
    if (prefix_sum(pc, 500ULL) != static_cast<u128>(13'343ULL)) {
        std::cerr << "Checkpoint failed: sum f(N), 1<=N<=500\n";
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

    const Precomp pc = build_precomp();

    if (options.run_checkpoints && !run_checkpoints(pc)) {
        return 1;
    }

    const u128 total = prefix_sum(pc, options.n);
    const u64 answer = static_cast<u64>(total % static_cast<u128>(kMod));
    std::cout << std::setw(8) << std::setfill('0') << answer << '\n';
    return 0;
}
