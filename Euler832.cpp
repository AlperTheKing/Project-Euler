#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kMod = 1'000'000'007ULL;
static constexpr u64 inv2 = (kMod + 1) / 2;

static u64 mod_add(u64 a, u64 b) {
    u64 s = a + b;
    if (s >= kMod) s -= kMod;
    return s;
}

static u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

static u64 range_sum_mod(u64 start, u64 len) {
    if (len == 0) return 0;
    u64 t1 = len % kMod;
    u64 t2 = (2 * (start % kMod) + (len - 1) % kMod) % kMod;
    return mod_mul(mod_mul(t1, t2), inv2);
}

static u64 full_block_sum(u64 a, u64 b, u64 c, u64 size) {
    u64 s = 0;
    s = mod_add(s, range_sum_mod(a, size));
    s = mod_add(s, range_sum_mod(b, size));
    s = mod_add(s, range_sum_mod(c, size));
    return s;
}

static u64 partial_block_sum(u64 len, u64 a, u64 b, u64 c, u64 size) {
    if (len == 0) return 0;
    if (len == size) return full_block_sum(a, b, c, size);
    if (size == 1) return (a + b + c) % kMod;

    const u64 s = size / 4;

    u64 ans = 0;

    u64 l1 = std::min(len, s);
    ans = mod_add(ans, partial_block_sum(l1, a, b, c, s));

    u64 rem = (len > s ? len - s : 0);
    u64 l2 = std::min(rem, s);
    ans = mod_add(ans, partial_block_sum(l2, a + s, b + 2 * s, c + 3 * s, s));

    rem = (len > 2 * s ? len - 2 * s : 0);
    u64 l3 = std::min(rem, s);
    ans = mod_add(ans, partial_block_sum(l3, a + 2 * s, b + 3 * s, c + s, s));

    rem = (len > 3 * s ? len - 3 * s : 0);
    u64 l4 = std::min(rem, s);
    ans = mod_add(ans, partial_block_sum(l4, a + 3 * s, b + s, c + 2 * s, s));

    return ans;
}

static u64 M_fast(u64 n) {
    if (n == 0) return 0;

    u64 a = 1, b = 2, c = 3;
    u64 size = 1;
    u64 left = n;
    u64 ans = 0;

    while (left > size) {
        ans = mod_add(ans, full_block_sum(a, b, c, size));
        left -= size;
        a *= 4;
        b *= 4;
        c *= 4;
        size *= 4;
    }

    ans = mod_add(ans, partial_block_sum(left, a, b, c, size));
    return ans;
}

static u64 M_bruteforce(int n) {
    std::set<int> paper;
    u64 sum = 0;

    for (int round = 0; round < n; ++round) {
        int a = 1;
        while (paper.count(a)) ++a;

        int b = 1;
        while (true) {
            int c = a ^ b;
            if (b > 0 && c > 0 && !paper.count(b) && !paper.count(c)) {
                break;
            }
            ++b;
        }

        int c = a ^ b;
        paper.insert(a);
        paper.insert(b);
        paper.insert(c);
    }

    for (int x : paper) sum += static_cast<u64>(x);
    return sum;
}

int main() {
    assert(M_bruteforce(10) == 642);
    assert(M_fast(10) == 642);
    assert(M_fast(1000) == 5432148);

    std::cout << M_fast(1'000'000'000'000'000'000ULL) << '\n';
    return 0;
}
