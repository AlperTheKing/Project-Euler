#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace {

using i128 = __int128_t;
using u64 = std::uint64_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 pow_u64(u64 base, int exp) {
    u64 result = 1;
    while (exp > 0) {
        if (exp & 1) {
            result *= base;
        }
        base *= base;
        exp >>= 1;
    }
    return result;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = static_cast<u64>((static_cast<i128>(result) * base) % kMod);
        }
        base = static_cast<u64>((static_cast<i128>(base) * base) % kMod);
        exp >>= 1ULL;
    }
    return result;
}

u64 G_mod(const int p) {
    const u64 A = pow_u64(17, p);
    const u64 B = pow_u64(19, p);
    const u64 C = pow_u64(23, p);

    const u64 inf = std::numeric_limits<u64>::max() / 4;
    std::vector<u64> dist(static_cast<std::size_t>(A), inf);

    using Node = std::pair<u64, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    dist[0] = 0;
    pq.push({0, 0});

    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[static_cast<std::size_t>(u)]) {
            continue;
        }

        int v = static_cast<int>((static_cast<u64>(u) + B) % A);
        u64 nd = d + B;
        if (nd < dist[static_cast<std::size_t>(v)]) {
            dist[static_cast<std::size_t>(v)] = nd;
            pq.push({nd, v});
        }

        v = static_cast<int>((static_cast<u64>(u) + C) % A);
        nd = d + C;
        if (nd < dist[static_cast<std::size_t>(v)]) {
            dist[static_cast<std::size_t>(v)] = nd;
            pq.push({nd, v});
        }
    }

    u64 sum_w_mod = 0;
    u64 sum_w2_mod = 0;
    for (const u64 w : dist) {
        const u64 wm = w % kMod;
        sum_w_mod += wm;
        if (sum_w_mod >= kMod) {
            sum_w_mod -= kMod;
        }
        sum_w2_mod = (sum_w2_mod + static_cast<u64>((static_cast<i128>(wm) * wm) % kMod)) % kMod;
    }

    const u64 inv2 = (kMod + 1) / 2;
    const u64 inv12 = mod_pow(12, kMod - 2);
    const u64 Am = A % kMod;
    const u64 invA = mod_pow(Am, kMod - 2);

    u64 genus = static_cast<u64>((static_cast<i128>(sum_w_mod) * invA) % kMod);
    genus = (genus + kMod - static_cast<u64>((static_cast<i128>((Am + kMod - 1) % kMod) * inv2) % kMod)) % kMod;

    u64 gaps = static_cast<u64>((static_cast<i128>(sum_w2_mod) * inv2) % kMod);
    gaps = static_cast<u64>((static_cast<i128>(gaps) * invA) % kMod);
    gaps = (gaps + kMod - static_cast<u64>((static_cast<i128>(sum_w_mod) * inv2) % kMod)) % kMod;
    const u64 A2m = static_cast<u64>((static_cast<i128>(Am) * Am) % kMod);
    gaps = (gaps + static_cast<u64>((static_cast<i128>((A2m + kMod - 1) % kMod) * inv12) % kMod)) % kMod;

    const u64 shift = (A % kMod + B % kMod + C % kMod) % kMod;
    const u64 tri = static_cast<u64>((static_cast<i128>(shift) * ((shift + kMod - 1) % kMod) % kMod) * inv2 % kMod);

    u64 ans = tri;
    ans = (ans + static_cast<u64>((static_cast<i128>(genus) * shift) % kMod)) % kMod;
    ans = (ans + gaps) % kMod;
    return ans;
}

}  // namespace

int main() {
    assert(G_mod(1) == 8'253);
    assert(G_mod(2) == 60'258'000);

    std::cout << G_mod(6) << '\n';
    return 0;
}
