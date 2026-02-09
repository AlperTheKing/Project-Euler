#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;

struct Options {
    u32 n = 2'000'000U;
    bool run_checkpoints = true;
};

struct Dir {
    int x = 0;
    int y = 0;
    u32 w = 0U;
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
        out = static_cast<u32>(std::stoul(tail));
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
    return options.n >= 3U;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const unsigned digit = static_cast<unsigned>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

int half_plane(const Dir& d) {
    return (d.y > 0 || (d.y == 0 && d.x > 0)) ? 0 : 1;
}

bool angle_less(const Dir& a, const Dir& b) {
    const int ha = half_plane(a);
    const int hb = half_plane(b);
    if (ha != hb) {
        return ha < hb;
    }
    const i64 cross = static_cast<i64>(a.x) * b.y - static_cast<i64>(a.y) * b.x;
    return cross > 0;
}

u128 solve(const u32 n) {
    std::vector<std::pair<int, int>> raw_dirs;
    raw_dirs.reserve(n);

    u32 x_mod = 1248U % 32323U;
    u32 y_mod = 8421U % 30103U;

    for (u32 i = 0U; i < n; ++i) {
        const int x = static_cast<int>(x_mod) - 16161;
        const int y = static_cast<int>(y_mod) - 15051;
        const int g = std::gcd(std::abs(x), std::abs(y));
        raw_dirs.emplace_back(x / g, y / g);

        x_mod = static_cast<u32>((static_cast<u64>(x_mod) * 1248ULL) % 32323ULL);
        y_mod = static_cast<u32>((static_cast<u64>(y_mod) * 8421ULL) % 30103ULL);
    }

    std::sort(raw_dirs.begin(), raw_dirs.end());

    std::vector<Dir> unique_lex;
    unique_lex.reserve(raw_dirs.size());
    for (const auto& p : raw_dirs) {
        if (!unique_lex.empty() && unique_lex.back().x == p.first && unique_lex.back().y == p.second) {
            ++unique_lex.back().w;
        } else {
            unique_lex.push_back(Dir{p.first, p.second, 1U});
        }
    }

    u128 W = 0;
    u128 S2 = 0;
    u128 S3 = 0;
    for (const Dir& d : unique_lex) {
        const u128 w = d.w;
        W += w;
        S2 += w * w;
        S3 += w * w * w;
    }

    const u128 total_distinct = (W * W * W - 3U * W * S2 + 2U * S3) / 6U;

    u128 boundary = 0;
    const auto lex_cmp = [](const Dir& a, const Dir& b) {
        if (a.x != b.x) {
            return a.x < b.x;
        }
        return a.y < b.y;
    };

    for (std::size_t i = 0; i < unique_lex.size(); ++i) {
        const Dir& d = unique_lex[i];
        const int ox = -d.x;
        const int oy = -d.y;

        Dir key{ox, oy, 0U};
        const auto it = std::lower_bound(unique_lex.begin(), unique_lex.end(), key, lex_cmp);
        if (it == unique_lex.end() || it->x != ox || it->y != oy) {
            continue;
        }

        if (d.x < ox || (d.x == ox && d.y < oy)) {
            const u128 wi = d.w;
            const u128 wj = it->w;
            boundary += wi * wj * (W - wi - wj);
        }
    }

    std::vector<Dir> by_angle = unique_lex;
    std::sort(by_angle.begin(), by_angle.end(), angle_less);
    const std::size_t m = by_angle.size();

    std::vector<u64> weight2(2 * m, 0ULL);
    for (std::size_t i = 0; i < m; ++i) {
        weight2[i] = by_angle[i].w;
        weight2[i + m] = by_angle[i].w;
    }

    std::vector<u64> pref_w(2 * m + 1, 0ULL);
    std::vector<u64> pref_w2(2 * m + 1, 0ULL);
    for (std::size_t i = 0; i < 2 * m; ++i) {
        pref_w[i + 1] = pref_w[i] + weight2[i];
        pref_w2[i + 1] = pref_w2[i] + weight2[i] * weight2[i];
    }

    std::vector<Dir> angle2(2 * m);
    for (std::size_t i = 0; i < m; ++i) {
        angle2[i] = by_angle[i];
        angle2[i + m] = by_angle[i];
    }

    u128 outside_open = 0;
    std::size_t j = 0;
    for (std::size_t i = 0; i < m; ++i) {
        if (j < i + 1) {
            j = i + 1;
        }
        const Dir& a = angle2[i];
        while (j < i + m) {
            const Dir& b = angle2[j];
            const i64 cross = static_cast<i64>(a.x) * b.y - static_cast<i64>(a.y) * b.x;
            if (cross > 0) {
                ++j;
            } else {
                break;
            }
        }

        const u128 S = static_cast<u128>(pref_w[j] - pref_w[i + 1]);
        const u128 SS = static_cast<u128>(pref_w2[j] - pref_w2[i + 1]);
        const u128 pair_weight = (S * S - SS) / 2U;
        outside_open += static_cast<u128>(by_angle[i].w) * pair_weight;
    }

    return total_distinct - outside_open - boundary;
}

bool run_checkpoints() {
    if (solve(8U) != 20U) {
        std::cerr << "Checkpoint failed: C(8)\n";
        return false;
    }
    if (solve(600U) != 8'950'634ULL) {
        std::cerr << "Checkpoint failed: C(600)\n";
        return false;
    }
    if (solve(40'000U) != 2'666'610'948'988ULL) {
        std::cerr << "Checkpoint failed: C(40000)\n";
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

    std::cout << to_string_u128(solve(options.n)) << '\n';
    return 0;
}
