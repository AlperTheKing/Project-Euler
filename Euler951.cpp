#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

struct Key {
    std::uint8_t r;
    std::uint8_t first;
    i64 a;
    i64 b;

    bool operator==(const Key& other) const {
        return r == other.r && first == other.first && a == other.a && b == other.b;
    }
};

struct KeyHash {
    std::size_t operator()(const Key& k) const {
        std::uint64_t h = 1469598103934665603ULL;
        auto mix = [&](std::uint64_t x) {
            h ^= x;
            h *= 1099511628211ULL;
        };
        mix(k.r);
        mix(k.first);
        mix(static_cast<std::uint64_t>(k.a));
        mix(static_cast<std::uint64_t>(k.b));
        return static_cast<std::size_t>(h);
    }
};

u64 F(int n) {
    const int m = 2 * n;

    std::unordered_map<Key, u64, KeyHash> cur;
    cur.reserve(2048);

    cur[{1, 1, 2, -1}] = 1ULL;
    cur[{0, 0, 2, -1}] = 1ULL;

    for (int len = 1; len < m; ++len) {
        std::unordered_map<Key, u64, KeyHash> nxt;
        nxt.reserve(cur.size() * 2 + 16);

        for (const auto& kv : cur) {
            const Key& key = kv.first;
            const u64 ways = kv.second;
            const int r = key.r;
            const int b = len - r;

            if (r < n) {
                const i64 ai = (key.first == 1) ? -(key.a + 2 * key.b) : (-2 * key.a);
                Key nk{static_cast<std::uint8_t>(r + 1), 1, ai, key.a};
                nxt[nk] += ways;
            }

            if (b < n) {
                const i64 ai = (key.first == 0) ? -(key.a + 2 * key.b) : (-2 * key.a);
                Key nk{static_cast<std::uint8_t>(r), 0, ai, key.a};
                nxt[nk] += ways;
            }
        }

        cur = std::move(nxt);
    }

    u64 ans = 0;
    for (const auto& kv : cur) {
        const Key& key = kv.first;
        if (key.r == n && key.a == 0) {
            ans += kv.second;
        }
    }

    return ans;
}

void run_validations() {
    assert(F(2) == 4ULL);
    assert(F(8) == 11'892ULL);
}

}  // namespace

int main() {
    run_validations();
    std::cout << F(26) << '\n';
    return 0;
}
