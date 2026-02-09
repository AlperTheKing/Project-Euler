#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

struct Options {
    u64 L = 1'000'000'000ULL;
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
        if (parse_u64_after_prefix(arg, "--L=", options.L)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

u128 tri_u128(const u64 n) {
    return static_cast<u128>(n) * static_cast<u128>(n + 1ULL) / 2U;
}

std::string to_string_u128(u128 v) {
    if (v == 0) {
        return "0";
    }
    std::string s;
    while (v > 0) {
        const int digit = static_cast<int>(v % 10U);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10U;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::vector<int> build_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            if (static_cast<int64_t>(i) * i <= n) {
                for (int x = i * i; x <= n; x += i) {
                    if (spf[static_cast<std::size_t>(x)] == 0) {
                        spf[static_cast<std::size_t>(x)] = i;
                    }
                }
            }
        }
    }
    spf[1] = 1;
    return spf;
}

std::vector<std::vector<std::pair<int, int>>> build_squarefree_mu_divisors(const int n,
                                                                            const std::vector<int>& spf) {
    std::vector<std::vector<std::pair<int, int>>> divs(static_cast<std::size_t>(n + 1));
    divs[1].push_back({1, 1});

    for (int x = 2; x <= n; ++x) {
        int t = x;
        std::vector<int> primes;
        while (t > 1) {
            const int p = spf[static_cast<std::size_t>(t)];
            primes.push_back(p);
            while (t % p == 0) {
                t /= p;
            }
        }

        std::vector<std::pair<int, int>> cur;
        cur.push_back({1, 1});
        for (const int p : primes) {
            const std::size_t before = cur.size();
            for (std::size_t i = 0; i < before; ++i) {
                cur.push_back({cur[i].first * p, -cur[i].second});
            }
        }
        divs[static_cast<std::size_t>(x)] = std::move(cur);
    }
    return divs;
}

i128 coprime_prefix_sum(const int r, const u64 x,
                        const std::vector<std::vector<std::pair<int, int>>>& sqf_divs) {
    if (x == 0ULL) {
        return 0;
    }
    i128 out = 0;
    for (const auto& [d, mu] : sqf_divs[static_cast<std::size_t>(r)]) {
        const u64 q = x / static_cast<u64>(d);
        const u128 tri = tri_u128(q);
        const i128 term = static_cast<i128>(static_cast<u128>(d) * tri);
        out += (mu > 0 ? term : -term);
    }
    return out;
}

u128 coprime_range_sum(const int r, const u64 l, const u64 rr,
                       const std::vector<std::vector<std::pair<int, int>>>& sqf_divs) {
    if (l > rr) {
        return 0;
    }
    const i128 hi = coprime_prefix_sum(r, rr, sqf_divs);
    const i128 lo = coprime_prefix_sum(r, l - 1ULL, sqf_divs);
    return static_cast<u128>(hi - lo);
}

u128 solve(const u64 L) {
    if (L < 2ULL) {
        return 0;
    }

    const int r_max = static_cast<int>(std::sqrt(static_cast<long double>(L)));
    const std::vector<int> spf = build_spf(r_max);
    const auto sqf_divs = build_squarefree_mu_divisors(r_max, spf);

    u128 total = 0;

    for (int r = 1; r <= r_max; ++r) {
        u64 s_lo = static_cast<u64>(r) + 1ULL;
        const u64 s_hi = std::min<u64>(2ULL * static_cast<u64>(r) - 1ULL, L / static_cast<u64>(r));
        if (s_lo > s_hi) {
            continue;
        }

        while (s_lo <= s_hi) {
            const u64 m = L / (static_cast<u64>(r) * s_lo);
            const u64 s_end =
                std::min<u64>(s_hi, L / (static_cast<u64>(r) * m));

            const u128 sum_s = coprime_range_sum(r, s_lo, s_end, sqf_divs);
            const u128 add = static_cast<u128>(r) * tri_u128(m) * sum_s;
            total += add;

            s_lo = s_end + 1ULL;
        }
    }

    return total;
}

bool run_checkpoints() {
    if (solve(15ULL) != static_cast<u128>(45ULL)) {
        std::cerr << "Checkpoint failed: F(15)\n";
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

    const u128 ans = solve(options.L);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
