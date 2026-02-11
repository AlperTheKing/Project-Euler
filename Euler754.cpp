#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kPhiMod = kMod - 1ULL;

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    u64 b = base % kMod;
    u64 e = exp;
    while (e > 0ULL) {
        if ((e & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<__uint128_t>(result) * b) % kMod);
        }
        b = static_cast<u64>((static_cast<__uint128_t>(b) * b) % kMod);
        e >>= 1ULL;
    }
    return result;
}

std::vector<std::int8_t> mobius_up_to(const int n) {
    std::vector<std::int8_t> mu(static_cast<std::size_t>(n + 1), 0);
    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(n + 1), 0U);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));

    mu[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (is_composite[static_cast<std::size_t>(i)] == 0U) {
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (const int p : primes) {
            const i64 v = static_cast<i64>(i) * static_cast<i64>(p);
            if (v > n) {
                break;
            }
            is_composite[static_cast<std::size_t>(v)] = 1U;
            if ((i % p) == 0) {
                mu[static_cast<std::size_t>(v)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(v)] = static_cast<std::int8_t>(-mu[static_cast<std::size_t>(i)]);
        }
    }

    return mu;
}

u64 G_mod(const int n) {
    const std::vector<std::int8_t> mu = mobius_up_to(n);

    struct Range {
        int l;
        int r;
        int q;
    };
    std::vector<Range> ranges;
    ranges.reserve(2U * static_cast<std::size_t>(std::sqrt(static_cast<long double>(n))) + 16U);

    std::vector<int> needed_k;
    needed_k.reserve(ranges.capacity() + 1U);
    needed_k.push_back(0);

    for (int l = 1; l <= n;) {
        const int q = n / l;
        const int r = n / q;
        ranges.push_back({l, r, q});
        needed_k.push_back(q - 1);
        l = r + 1;
    }

    std::sort(needed_k.begin(), needed_k.end());
    needed_k.erase(std::unique(needed_k.begin(), needed_k.end()), needed_k.end());

    std::unordered_map<int, int> idx;
    idx.reserve(needed_k.size() * 2U);
    for (int i = 0; i < static_cast<int>(needed_k.size()); ++i) {
        idx[needed_k[static_cast<std::size_t>(i)]] = i;
    }

    std::vector<u64> superfac_value(needed_k.size(), 0ULL);
    superfac_value[static_cast<std::size_t>(idx[0])] = 1ULL;

    const int max_k = needed_k.back();
    u64 fact = 1ULL;
    u64 superfac = 1ULL;
    std::size_t ptr = 0;
    if (needed_k[ptr] == 0) {
        ++ptr;
    }

    for (int t = 1; t <= max_k; ++t) {
        fact = static_cast<u64>((static_cast<__uint128_t>(fact) * static_cast<u64>(t)) % kMod);
        superfac = static_cast<u64>((static_cast<__uint128_t>(superfac) * fact) % kMod);
        while (ptr < needed_k.size() && needed_k[ptr] == t) {
            superfac_value[static_cast<std::size_t>(ptr)] = superfac;
            ++ptr;
        }
    }

    u64 answer = 1ULL;

    for (const Range& rg : ranges) {
        const u64 q = static_cast<u64>(rg.q);
        const u64 f = ((q * (q - 1ULL)) / 2ULL) % kPhiMod;

        u64 prod_pos = 1ULL;
        u64 prod_neg = 1ULL;
        i64 sum_mu = 0;

        for (int d = rg.l; d <= rg.r; ++d) {
            const std::int8_t m = mu[static_cast<std::size_t>(d)];
            if (m == 0) {
                continue;
            }
            sum_mu += static_cast<i64>(m);
            if (m > 0) {
                prod_pos = static_cast<u64>((static_cast<__uint128_t>(prod_pos) * static_cast<u64>(d)) % kMod);
            } else {
                prod_neg = static_cast<u64>((static_cast<__uint128_t>(prod_neg) * static_cast<u64>(d)) % kMod);
            }
        }

        if (f != 0ULL) {
            answer = static_cast<u64>((static_cast<__uint128_t>(answer) * mod_pow(prod_pos, f)) % kMod);
            answer = static_cast<u64>(
                (static_cast<__uint128_t>(answer) * mod_pow(prod_neg, (kPhiMod - f) % kPhiMod)) % kMod);
        }

        i64 exp_sf = sum_mu % static_cast<i64>(kPhiMod);
        if (exp_sf < 0) {
            exp_sf += static_cast<i64>(kPhiMod);
        }
        const int key = rg.q - 1;
        const int k_index = idx[key];
        const u64 sf = superfac_value[static_cast<std::size_t>(k_index)];
        answer = static_cast<u64>((static_cast<__uint128_t>(answer) * mod_pow(sf, static_cast<u64>(exp_sf))) % kMod);
    }

    return answer;
}

}  // namespace

int main() {
    assert(G_mod(10) == 331'358'692ULL);
    std::cout << G_mod(100'000'000) << '\n';
    return 0;
}
