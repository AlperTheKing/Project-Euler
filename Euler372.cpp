#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 isqrt_u128(const u128 n) {
    if (n == 0U) {
        return 0U;
    }
    long double approx = std::sqrt(static_cast<long double>(n));
    u64 r = static_cast<u64>(approx);
    while (static_cast<u128>(r + 1U) * static_cast<u128>(r + 1U) <= n) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > n) {
        --r;
    }
    return r;
}

u64 floor_div_by_sqrt(const u64 t, const int d) {
    const u128 v = (static_cast<u128>(t) * static_cast<u128>(t)) / static_cast<u128>(d);
    return isqrt_u128(v);
}

struct BeattyKey {
    int d;
    int p;
    int q;
    u64 n;

    bool operator==(const BeattyKey& other) const {
        return d == other.d && p == other.p && q == other.q && n == other.n;
    }
};

struct BeattyKeyHash {
    std::size_t operator()(const BeattyKey& k) const {
        std::size_t h = static_cast<std::size_t>(k.d);
        h = h * 1315423911u + static_cast<std::size_t>(k.p);
        h = h * 1315423911u + static_cast<std::size_t>(k.q);
        h = h * 1315423911u + static_cast<std::size_t>(k.n ^ (k.n >> 32U));
        return h;
    }
};

class BeattySumSqrt {
public:
    BeattySumSqrt() = default;

    u128 sum_floor(const int d, const u64 n) {
        return sum_floor_impl(d, n, 0, 1);
    }

private:
    std::unordered_map<BeattyKey, u128, BeattyKeyHash> memo_;

    static bool ge_asqrt_d(const u64 a, const int d, const i64 rhs) {
        if (rhs <= 0) {
            return true;
        }
        const u128 left = static_cast<u128>(a) * static_cast<u128>(a) * static_cast<u128>(d);
        const u128 right = static_cast<u128>(rhs) * static_cast<u128>(rhs);
        return left >= right;
    }

    static i64 floor_linear_surd(const u64 a, const int d, const i64 b, const int c) {
        // floor((a*sqrt(d) + b) / c), d is non-square.
        long double guess_ld =
            (static_cast<long double>(a) * std::sqrt(static_cast<long double>(d)) + static_cast<long double>(b)) /
            static_cast<long double>(c);
        i64 g = static_cast<i64>(std::floor(guess_ld));

        while (ge_asqrt_d(a, d, static_cast<i64>((g + 1) * static_cast<i64>(c) - b))) {
            ++g;
        }
        while (!ge_asqrt_d(a, d, static_cast<i64>(g * static_cast<i64>(c) - b))) {
            --g;
        }
        return g;
    }

    u128 sum_floor_impl(const int d, const u64 n, const int p, const int q) {
        if (n == 0U) {
            return 0U;
        }

        const int s = static_cast<int>(std::sqrt(static_cast<long double>(d)));
        if (s * s == d) {
            return static_cast<u128>(s) * static_cast<u128>(n) * static_cast<u128>(n + 1U) / 2U;
        }

        const BeattyKey key{d, p, q, n};
        const auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        const int a = (s + p) / q;
        const int p1 = a * q - p;
        const int q1 = (d - p1 * p1) / q;

        const i64 m_signed = floor_linear_surd(n, d, -static_cast<i64>(n) * p1, q);
        const u64 m = static_cast<u64>(m_signed);

        const u128 tri = static_cast<u128>(n) * static_cast<u128>(n + 1U) / 2U;
        const u128 result =
            static_cast<u128>(a) * tri + static_cast<u128>(n) * static_cast<u128>(m) - sum_floor_impl(d, m, p1, q1);

        memo_.emplace(key, result);
        return result;
    }
};

u128 sum_ceil_sqrt(BeattySumSqrt& beatty, const int d, const u64 l, const u64 r) {
    if (l > r) {
        return 0U;
    }

    const int s = static_cast<int>(std::sqrt(static_cast<long double>(d)));
    if (s * s == d) {
        const u128 cnt = static_cast<u128>(r - l + 1U);
        return static_cast<u128>(s) * static_cast<u128>(l + r) * cnt / 2U;
    }

    const u128 floors = beatty.sum_floor(d, r) - beatty.sum_floor(d, l - 1U);
    return floors + static_cast<u128>(r - l + 1U);
}

u128 count_rays(const u64 m, const u64 n) {
    const u64 l = m + 1U;
    const u64 k_max = static_cast<u64>((static_cast<u128>(n) * static_cast<u128>(n)) /
                                       (static_cast<u128>(l) * static_cast<u128>(l)));

    BeattySumSqrt beatty;
    u128 total = 0U;

    for (u64 k = 1U; k <= k_max; k += 2U) {
        const u64 u = floor_div_by_sqrt(n, static_cast<int>(k));
        if (u < l) {
            continue;
        }

        const u64 v = floor_div_by_sqrt(n + 1U, static_cast<int>(k + 1U));
        const u64 x1 = std::min(u, v);
        if (x1 >= l) {
            total += sum_ceil_sqrt(beatty, static_cast<int>(k + 1U), l, x1) -
                     sum_ceil_sqrt(beatty, static_cast<int>(k), l, x1);
        }

        const u64 x2 = std::max(l, v + 1U);
        if (u >= x2) {
            total += static_cast<u128>(u - x2 + 1U) * static_cast<u128>(n + 1U) -
                     sum_ceil_sqrt(beatty, static_cast<int>(k), x2, u);
        }
    }

    return total;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }
    std::string out;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool run_checkpoints() {
    if (count_rays(0U, 100U) != 3019U) {
        std::cerr << "Checkpoint failed: R(0,100)\n";
        return false;
    }
    if (count_rays(100U, 10000U) != 29750422U) {
        std::cerr << "Checkpoint failed: R(100,10000)\n";
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

    const u128 answer = count_rays(2000000ULL, 1000000000ULL);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
