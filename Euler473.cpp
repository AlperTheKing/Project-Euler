#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <functional>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    i64 limit = 10'000'000'000LL;
    bool run_checkpoints = true;
};

struct LowState {
    i64 b = 0;
    i64 a = 0;
    std::uint8_t last = 0;
};

struct HighState {
    i64 b = 0;
    i64 a = 0;
    std::uint8_t first = 0;
};

struct Bucket {
    std::vector<i64> values;
    std::vector<i64> pref;  // pref[0]=0, pref[i+1]=sum values[0..i]
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<i64>(std::stoll(tail));
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.limit <= 0) {
        std::cerr << "--limit must be positive.\n";
        return false;
    }
    return true;
}

std::array<i64, 100> build_fib() {
    std::array<i64, 100> f{};
    f[0] = 0;
    f[1] = 1;
    for (int i = 2; i < 100; ++i) {
        f[i] = f[i - 1] + f[i - 2];
    }
    return f;
}

void build_pair_coefficients(const std::array<i64, 100>& fib,
                             std::vector<i64>& coeff_a,
                             std::vector<i64>& coeff_b,
                             const int max_index) {
    coeff_a.assign(static_cast<std::size_t>(max_index + 1), 0);
    coeff_b.assign(static_cast<std::size_t>(max_index + 1), 0);

    for (int i = 0; i <= max_index; ++i) {
        i64 pos_a = 0;
        i64 pos_b = 0;
        if (i == 0) {
            pos_a = 1;
            pos_b = 0;
        } else {
            pos_a = fib[static_cast<std::size_t>(i - 1)];
            pos_b = fib[static_cast<std::size_t>(i)];
        }

        const int m = i + 1;
        const i64 sign = (m % 2 == 0) ? 1 : -1;  // (-1)^m
        const i64 neg_a = sign * fib[static_cast<std::size_t>(m + 1)];
        const i64 neg_b = -sign * fib[static_cast<std::size_t>(m)];

        coeff_a[static_cast<std::size_t>(i)] = pos_a + neg_a;
        coeff_b[static_cast<std::size_t>(i)] = pos_b + neg_b;
    }
}

std::pair<i64, i64> min_max_possible_a_for_length(const int n, const std::vector<i64>& coeff_a) {
    constexpr i64 kInf = (1LL << 62);
    std::vector<std::array<i64, 2>> min_dp(static_cast<std::size_t>(n + 1), {kInf, kInf});
    std::vector<std::array<i64, 2>> max_dp(static_cast<std::size_t>(n + 1), {-kInf, -kInf});
    min_dp[static_cast<std::size_t>(n)] = {0, 0};
    max_dp[static_cast<std::size_t>(n)] = {0, 0};

    for (int idx = n - 1; idx >= 0; --idx) {
        for (int prev = 0; prev <= 1; ++prev) {
            i64 best_min = kInf;
            i64 best_max = -kInf;

            auto relax = [&](const int bit) {
                if (prev == 1 && bit == 1) {
                    return;
                }
                const i64 add = (bit == 1) ? coeff_a[static_cast<std::size_t>(idx)] : 0;
                const i64 child_min = min_dp[static_cast<std::size_t>(idx + 1)][bit];
                const i64 child_max = max_dp[static_cast<std::size_t>(idx + 1)][bit];
                if (child_min >= kInf / 2 || child_max <= -kInf / 2) {
                    return;
                }
                best_min = std::min(best_min, add + child_min);
                best_max = std::max(best_max, add + child_max);
            };

            if (idx == 0) {
                relax(0);  // d_0 = 0 (otherwise d_0 and d_{-1} would be consecutive 1's in palindrome)
            } else if (idx == n - 1) {
                relax(1);  // highest digit must be 1 (no leading zero)
            } else {
                relax(0);
                relax(1);
            }

            min_dp[static_cast<std::size_t>(idx)][prev] = best_min;
            max_dp[static_cast<std::size_t>(idx)][prev] = best_max;
        }
    }

    return {min_dp[0][0], max_dp[0][0]};
}

void enumerate_low_states(const int end,
                          const std::vector<i64>& coeff_a,
                          const std::vector<i64>& coeff_b,
                          std::vector<LowState>& out) {
    out.clear();
    out.reserve(1U << 16);

    std::function<void(int, int, i64, i64, std::uint8_t)> dfs =
        [&](const int idx, const int prev, const i64 sum_b, const i64 sum_a, const std::uint8_t last) {
            if (idx > end) {
                out.push_back(LowState{sum_b, sum_a, last});
                return;
            }
            if (idx == 0) {
                dfs(1, 0, sum_b, sum_a, 0);
                return;
            }

            // bit = 0
            dfs(idx + 1, 0, sum_b, sum_a, 0);

            // bit = 1
            if (prev == 0) {
                dfs(idx + 1, 1, sum_b + coeff_b[static_cast<std::size_t>(idx)],
                    sum_a + coeff_a[static_cast<std::size_t>(idx)], 1);
            }
        };

    if (end < 0) {
        out.push_back(LowState{0, 0, 0});
        return;
    }
    dfs(0, 0, 0, 0, 0);
}

void enumerate_high_states(const int start,
                           const int end,
                           const std::vector<i64>& coeff_a,
                           const std::vector<i64>& coeff_b,
                           std::vector<HighState>& out) {
    out.clear();
    out.reserve(1U << 16);

    std::function<void(int, int, i64, i64, std::uint8_t)> dfs =
        [&](const int idx, const int prev, const i64 sum_b, const i64 sum_a, const std::uint8_t first) {
            if (idx > end) {
                out.push_back(HighState{sum_b, sum_a, first});
                return;
            }

            auto first_after = [&](const int bit) -> std::uint8_t {
                return (idx == start) ? static_cast<std::uint8_t>(bit) : first;
            };

            if (idx == end) {
                const int bit = 1;  // highest digit fixed to 1
                if (prev == 1) {
                    return;
                }
                dfs(idx + 1, 1, sum_b + coeff_b[static_cast<std::size_t>(idx)],
                    sum_a + coeff_a[static_cast<std::size_t>(idx)], first_after(1));
                return;
            }

            // bit = 0
            dfs(idx + 1, 0, sum_b, sum_a, first_after(0));

            // bit = 1
            if (prev == 0) {
                dfs(idx + 1, 1, sum_b + coeff_b[static_cast<std::size_t>(idx)],
                    sum_a + coeff_a[static_cast<std::size_t>(idx)], first_after(1));
            }
        };

    if (start > end) {
        out.push_back(HighState{0, 0, 0});
        return;
    }
    dfs(start, 0, 0, 0, 0);
}

u128 solve_sum(const i64 limit) {
    const int max_index = 64;
    const auto fib = build_fib();
    std::vector<i64> coeff_a;
    std::vector<i64> coeff_b;
    build_pair_coefficients(fib, coeff_a, coeff_b, max_index);

    u128 total = 1;  // representation "1"

    std::vector<LowState> low_states;
    std::vector<HighState> high_states;

    for (int n = 2; n <= 64; ++n) {
        const auto [min_a, max_a] = min_max_possible_a_for_length(n, coeff_a);
        if (max_a < 1 || min_a > limit) {
            continue;
        }

        const int mid = n / 2;
        enumerate_low_states(mid - 1, coeff_a, coeff_b, low_states);
        enumerate_high_states(mid, n - 1, coeff_a, coeff_b, high_states);

        std::array<std::unordered_map<i64, Bucket>, 2> high_maps;
        for (int fb = 0; fb <= 1; ++fb) {
            high_maps[fb].reserve(high_states.size() / 2U + 8U);
        }

        for (const HighState& st : high_states) {
            high_maps[st.first][st.b].values.push_back(st.a);
        }
        for (int fb = 0; fb <= 1; ++fb) {
            for (auto& kv : high_maps[fb]) {
                auto& vals = kv.second.values;
                std::sort(vals.begin(), vals.end());
                auto& pref = kv.second.pref;
                pref.assign(vals.size() + 1U, 0);
                for (std::size_t i = 0; i < vals.size(); ++i) {
                    pref[i + 1U] = pref[i] + vals[i];
                }
            }
        }

        for (const LowState& st : low_states) {
            for (int fb = 0; fb <= 1; ++fb) {
                if (st.last == 1U && fb == 1) {
                    continue;
                }
                const i64 target_b = -st.b;
                auto it = high_maps[fb].find(target_b);
                if (it == high_maps[fb].end()) {
                    continue;
                }

                const auto& vals = it->second.values;
                const auto& pref = it->second.pref;
                const i64 low_needed = 1 - st.a;
                const i64 high_needed = limit - st.a;

                const auto lo_it = std::lower_bound(vals.begin(), vals.end(), low_needed);
                const auto hi_it = std::upper_bound(vals.begin(), vals.end(), high_needed);
                if (lo_it == hi_it) {
                    continue;
                }

                const std::size_t lo = static_cast<std::size_t>(lo_it - vals.begin());
                const std::size_t hi = static_cast<std::size_t>(hi_it - vals.begin());
                const u128 count = static_cast<u128>(hi - lo);
                const i64 sum_high = pref[hi] - pref[lo];
                total += count * static_cast<u128>(st.a) + static_cast<u128>(sum_high);
            }
        }
    }

    return total;
}

bool run_checkpoints() {
    if (solve_sum(1000LL) != static_cast<u128>(4345ULL)) {
        std::cerr << "Checkpoint failed: sum <= 1000\n";
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

    const u128 answer = solve_sum(options.limit);
    std::cout << static_cast<u64>(answer) << '\n';
    return 0;
}
