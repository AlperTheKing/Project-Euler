#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;

struct Node {
    i64 sum = 0;
    i64 pref = 0;
    i64 suff = 0;
    i64 best = 0;
};

Node merge_nodes(const Node& a, const Node& b) {
    Node c;
    c.sum = a.sum + b.sum;
    c.pref = std::max(a.pref, a.sum + b.pref);
    c.suff = std::max(b.suff, b.sum + a.suff);
    c.best = std::max({a.best, b.best, a.suff + b.pref});
    return c;
}

class MaxSubarrayTracker {
public:
    MaxSubarrayTracker(int n, int block_size) : n_(n), bsz_(block_size) {
        block_count_ = (n_ + bsz_ - 1) / bsz_;
        arr_.assign(static_cast<std::size_t>(n_), 0);
        size_ = 1;
        while (size_ < block_count_) size_ <<= 1;
        seg_.assign(static_cast<std::size_t>(2 * size_), Node{});
    }

    void add(int idx, i64 delta) {
        arr_[static_cast<std::size_t>(idx)] += delta;
        const int b = idx / bsz_;
        const Node leaf = recompute_block(b);
        int pos = size_ + b;
        seg_[static_cast<std::size_t>(pos)] = leaf;
        pos >>= 1;
        while (pos > 0) {
            seg_[static_cast<std::size_t>(pos)] =
                merge_nodes(seg_[static_cast<std::size_t>(pos << 1)],
                            seg_[static_cast<std::size_t>((pos << 1) | 1)]);
            pos >>= 1;
        }
    }

    i64 max_subarray() const { return seg_[1].best; }

private:
    Node recompute_block(int b) const {
        const int l = b * bsz_;
        const int r = std::min(n_, l + bsz_);

        Node x;
        i64 cur = 0;
        i64 run_pref = 0;
        for (int i = l; i < r; ++i) {
            const i64 v = arr_[static_cast<std::size_t>(i)];
            x.sum += v;
            run_pref += v;
            if (run_pref > x.pref) x.pref = run_pref;

            cur = std::max<i64>(0, cur + v);
            if (cur > x.best) x.best = cur;
        }

        i64 run_suff = 0;
        for (int i = r - 1; i >= l; --i) {
            run_suff += arr_[static_cast<std::size_t>(i)];
            if (run_suff > x.suff) x.suff = run_suff;
        }
        return x;
    }

    int n_;
    int bsz_;
    int block_count_;
    int size_;
    std::vector<i64> arr_;
    std::vector<Node> seg_;
};

i128 solve_range(int n, int l_from_exclusive, int l_to_inclusive) {
    const int block_size = 64;
    MaxSubarrayTracker tracker(n, block_size);

    const int need = 2 * l_to_inclusive;
    std::vector<int> t(static_cast<std::size_t>(need), 0);
    t[0] = 0;
    t[1] = 0;
    t[2] = 1 % n;
    for (int k = 3; k < need; ++k) {
        const int v = t[static_cast<std::size_t>(k - 1)] + t[static_cast<std::size_t>(k - 2)] +
                      t[static_cast<std::size_t>(k - 3)];
        t[static_cast<std::size_t>(k)] = v % n;
    }

    i128 answer = 0;
    for (int i = 1; i <= l_to_inclusive; ++i) {
        const int idx = t[static_cast<std::size_t>(2 * i - 2)];
        const i64 delta = 2LL * static_cast<i64>(t[static_cast<std::size_t>(2 * i - 1)]) -
                          static_cast<i64>(n) + 1LL;
        tracker.add(idx, delta);
        if (i > l_from_exclusive) answer += static_cast<i128>(tracker.max_subarray());
    }
    return answer;
}

i128 S(int n, int l) { return solve_range(n, 0, l); }

std::string to_string_i128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

}  // namespace

int main() {
    assert(S(5, 6) == 32);
    assert(S(5, 100) == 2416);
    assert(S(14, 100) == 3881);
    assert(S(107, 1000) == 1618572);

    const i128 ans = solve_range(10'000'003, 10'000'000, 10'200'000);
    std::cout << to_string_i128(ans) << "\n";
    return 0;
}
