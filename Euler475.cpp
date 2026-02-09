#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using i64 = std::int64_t;

constexpr u64 kMod = 1'000'000'007ULL;

struct Options {
    int musicians = 600;  // 12n
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--musicians=", options.musicians)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.musicians <= 0 || options.musicians % 12 != 0) {
        std::cerr << "--musicians must be a positive multiple of 12.\n";
        return false;
    }
    return true;
}

u64 mod_pow(u64 base, i64 exp) {
    u64 out = 1ULL;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1LL) {
            out = static_cast<u64>((__uint128_t)out * base % kMod);
        }
        base = static_cast<u64>((__uint128_t)base * base % kMod);
        exp >>= 1LL;
    }
    return out;
}

u32 pack_state(const int c0, const int c1, const int c2, const int c3) {
    return static_cast<u32>((c0 << 24) | (c1 << 16) | (c2 << 8) | c3);
}

void unpack_state(const u32 key, int& c0, int& c1, int& c2, int& c3) {
    c0 = static_cast<int>((key >> 24) & 0xFFU);
    c1 = static_cast<int>((key >> 16) & 0xFFU);
    c2 = static_cast<int>((key >> 8) & 0xFFU);
    c3 = static_cast<int>(key & 0xFFU);
}

u64 solve(const int musicians) {
    const int n = musicians / 12;
    const int m = 3 * n;   // quartets on day 1
    const int t = 4 * n;   // trios on day 2

    std::vector<std::array<u64, 4>> comb(static_cast<std::size_t>(m + 1));
    for (int i = 0; i <= m; ++i) {
        comb[static_cast<std::size_t>(i)][0] = 1ULL;
        if (i >= 1) {
            comb[static_cast<std::size_t>(i)][1] = static_cast<u64>(i);
        }
        if (i >= 2) {
            comb[static_cast<std::size_t>(i)][2] = static_cast<u64>(i) * (i - 1) / 2ULL;
        }
        if (i >= 3) {
            comb[static_cast<std::size_t>(i)][3] =
                static_cast<u64>(i) * (i - 1) * (i - 2) / 6ULL;
        }
    }

    std::vector<std::array<int, 4>> picks;
    picks.reserve(20);
    for (int x0 = 0; x0 <= 3; ++x0) {
        for (int x1 = 0; x1 + x0 <= 3; ++x1) {
            for (int x2 = 0; x2 + x1 + x0 <= 3; ++x2) {
                const int x3 = 3 - x0 - x1 - x2;
                picks.push_back({x0, x1, x2, x3});
            }
        }
    }

    std::unordered_map<u32, u64> cur;
    std::unordered_map<u32, u64> nxt;
    cur.reserve(200'000U);
    nxt.reserve(200'000U);
    cur[pack_state(m, 0, 0, 0)] = 1ULL;

    for (int step = 0; step < t; ++step) {
        nxt.clear();
        for (const auto& kv : cur) {
            int c0 = 0;
            int c1 = 0;
            int c2 = 0;
            int c3 = 0;
            unpack_state(kv.first, c0, c1, c2, c3);
            const u64 ways_so_far = kv.second;

            for (const auto& pick : picks) {
                const int x0 = pick[0];
                const int x1 = pick[1];
                const int x2 = pick[2];
                const int x3 = pick[3];
                if (x0 > c0 || x1 > c1 || x2 > c2 || x3 > c3) {
                    continue;
                }

                const int nc0 = c0 - x0;
                const int nc1 = c1 - x1 + x0;
                const int nc2 = c2 - x2 + x1;
                const int nc3 = c3 - x3 + x2;

                u64 ways_pick = comb[static_cast<std::size_t>(c0)][static_cast<std::size_t>(x0)];
                ways_pick = static_cast<u64>((__uint128_t)ways_pick *
                                             comb[static_cast<std::size_t>(c1)]
                                                 [static_cast<std::size_t>(x1)] %
                                             kMod);
                ways_pick = static_cast<u64>((__uint128_t)ways_pick *
                                             comb[static_cast<std::size_t>(c2)]
                                                 [static_cast<std::size_t>(x2)] %
                                             kMod);
                ways_pick = static_cast<u64>((__uint128_t)ways_pick *
                                             comb[static_cast<std::size_t>(c3)]
                                                 [static_cast<std::size_t>(x3)] %
                                             kMod);

                const u64 add =
                    static_cast<u64>((__uint128_t)ways_so_far * ways_pick % kMod);
                const u32 nkey = pack_state(nc0, nc1, nc2, nc3);
                auto it = nxt.find(nkey);
                if (it == nxt.end()) {
                    nxt.emplace(nkey, add);
                } else {
                    u64 v = it->second + add;
                    if (v >= kMod) {
                        v -= kMod;
                    }
                    it->second = v;
                }
            }
        }
        cur.swap(nxt);
    }

    const u32 final_key = pack_state(0, 0, 0, 0);
    const u64 matrix_count = cur[final_key];

    const u64 assign_quartet_members = mod_pow(24ULL, m);
    u64 fact_t = 1ULL;
    for (int i = 2; i <= t; ++i) {
        fact_t = static_cast<u64>((__uint128_t)fact_t * static_cast<u64>(i) % kMod);
    }
    const u64 inv_fact_t = mod_pow(fact_t, static_cast<i64>(kMod - 2ULL));

    u64 out = static_cast<u64>((__uint128_t)matrix_count * assign_quartet_members % kMod);
    out = static_cast<u64>((__uint128_t)out * inv_fact_t % kMod);
    return out;
}

bool run_checkpoints() {
    if (solve(12) != 576ULL) {
        std::cerr << "Checkpoint failed: f(12)\n";
        return false;
    }
    if (solve(24) != 509'089'824ULL) {
        std::cerr << "Checkpoint failed: f(24)\n";
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
        return 1;
    }

    std::cout << solve(options.musicians) << '\n';
    return 0;
}
