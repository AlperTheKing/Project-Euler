#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <queue>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kMod = 1'234'567'891ULL;

struct Node {
    long double logv;
    u64 modv;
    int id;
};

struct MinCmp {
    bool operator()(const Node& a, const Node& b) const {
        if (a.logv != b.logv) return a.logv > b.logv;
        return a.id > b.id;
    }
};

static u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 res = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = static_cast<u64>((static_cast<u128>(res) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return res;
}

static u64 S(u64 n, u64 m) {
    const int len = static_cast<int>(n - 1);

    std::priority_queue<Node, std::vector<Node>, MinCmp> pq;
    std::vector<Node> cur(len);
    long double max_log = 0.0L;

    for (int i = 0; i < len; ++i) {
        const u64 v = static_cast<u64>(i + 2);
        cur[i] = {std::log(static_cast<long double>(v)), v % kMod, i};
        pq.push(cur[i]);
        if (cur[i].logv > max_log) {
            max_log = cur[i].logv;
        }
    }

    u64 steps = 0;
    while (steps < m) {
        Node x = pq.top();
        if (2.0L * x.logv > max_log) {
            break;
        }
        pq.pop();
        x.logv *= 2.0L;
        x.modv = static_cast<u64>((static_cast<u128>(x.modv) * x.modv) % kMod);
        pq.push(x);
        if (x.logv > max_log) {
            max_log = x.logv;
        }
        ++steps;
    }

    std::vector<Node> arr;
    arr.reserve(len);
    while (!pq.empty()) {
        arr.push_back(pq.top());
        pq.pop();
    }

    if (steps == m) {
        u64 ans = 0;
        for (const auto& x : arr) {
            ans += x.modv;
            if (ans >= kMod) {
                ans -= kMod;
            }
        }
        return ans;
    }

    std::sort(arr.begin(), arr.end(), [](const Node& a, const Node& b) {
        if (a.logv != b.logv) return a.logv < b.logv;
        return a.id < b.id;
    });

    const u64 rem = m - steps;
    const u64 q = rem / static_cast<u64>(len);
    const u64 r = rem % static_cast<u64>(len);

    const u64 exp_q = mod_pow(2, q, kMod - 1);
    const u64 exp_q1 = static_cast<u64>((static_cast<u128>(exp_q) * 2ULL) % (kMod - 1));

    u64 ans = 0;
    for (int i = 0; i < len; ++i) {
        const u64 e = (static_cast<u64>(i) < r) ? exp_q1 : exp_q;
        const u64 add = mod_pow(arr[i].modv, e, kMod);
        ans += add;
        if (ans >= kMod) {
            ans -= kMod;
        }
    }

    return ans;
}

int main() {
    assert(S(5, 3) == 34);
    assert(S(10, 100) == 845339386ULL);
    std::cout << S(10'000, 10'000'000'000'000'000ULL) << '\n';
    return 0;
}
