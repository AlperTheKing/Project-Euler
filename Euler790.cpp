#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;

static constexpr i64 MOD = 50'515'093;

struct Rect {
    int x1;
    int x2;
    int y1;
    int y2;
};

struct Event {
    int x;
    int y1;
    int y2;
    int delta;
    bool operator<(const Event& other) const {
        return x < other.x;
    }
};

struct SegTree {
    int n;
    std::vector<std::array<i64, 12>> sum;
    std::vector<int> lazy;

    explicit SegTree(const std::vector<int>& ys) {
        n = static_cast<int>(ys.size()) - 1;
        sum.assign(4 * n + 4, {});
        lazy.assign(4 * n + 4, 0);
        build(1, 0, n - 1, ys);
    }

    static void rotate_arr(std::array<i64, 12>& a, int shift) {
        shift %= 12;
        if (shift < 0) shift += 12;
        if (shift == 0) return;
        std::array<i64, 12> b{};
        for (int i = 0; i < 12; ++i) {
            b[(i + shift) % 12] = a[i];
        }
        a = b;
    }

    void apply(int idx, int shift) {
        rotate_arr(sum[idx], shift);
        lazy[idx] += shift;
        lazy[idx] %= 12;
        if (lazy[idx] < 0) lazy[idx] += 12;
    }

    void pull(int idx) {
        for (int r = 0; r < 12; ++r) {
            sum[idx][r] = sum[idx * 2][r] + sum[idx * 2 + 1][r];
        }
    }

    void push(int idx) {
        if (lazy[idx] != 0) {
            apply(idx * 2, lazy[idx]);
            apply(idx * 2 + 1, lazy[idx]);
            lazy[idx] = 0;
        }
    }

    void build(int idx, int l, int r, const std::vector<int>& ys) {
        if (l == r) {
            const i64 len = static_cast<i64>(ys[l + 1]) - static_cast<i64>(ys[l]);
            sum[idx][0] = len;
            return;
        }
        const int m = (l + r) / 2;
        build(idx * 2, l, m, ys);
        build(idx * 2 + 1, m + 1, r, ys);
        pull(idx);
    }

    void update(int idx, int l, int r, int ql, int qr, int delta) {
        if (ql > r || qr < l) return;
        if (ql <= l && r <= qr) {
            apply(idx, delta);
            return;
        }
        push(idx);
        const int m = (l + r) / 2;
        update(idx * 2, l, m, ql, qr, delta);
        update(idx * 2 + 1, m + 1, r, ql, qr, delta);
        pull(idx);
    }

    void update(int ql, int qr, int delta) {
        if (ql > qr) return;
        update(1, 0, n - 1, ql, qr, delta);
    }

    const std::array<i64, 12>& all() const {
        return sum[1];
    }
};

static std::vector<Rect> generate_rects(int t) {
    std::vector<Rect> out;
    out.reserve(t);

    i64 s = 290'797;
    auto next_s = [&]() {
        s = (s * s) % MOD;
        return static_cast<int>(s);
    };

    std::vector<int> vals;
    vals.reserve(4 * t);
    vals.push_back(static_cast<int>(s));
    for (int i = 1; i < 4 * t; ++i) {
        vals.push_back(next_s());
    }

    for (int i = 0; i < t; ++i) {
        const int a = vals[4 * i + 0];
        const int b = vals[4 * i + 1];
        const int c = vals[4 * i + 2];
        const int d = vals[4 * i + 3];

        const int x1 = std::min(a, b);
        const int x2 = std::max(a, b);
        const int y1 = std::min(c, d);
        const int y2 = std::max(c, d);
        out.push_back({x1, x2, y1, y2});
    }

    return out;
}

static u64 solve(int t) {
    if (t == 0) {
        return static_cast<u64>(MOD) * static_cast<u64>(MOD) * 12ULL;
    }

    const auto rects = generate_rects(t);

    std::vector<int> ys;
    ys.reserve(2 * t + 2);
    ys.push_back(0);
    ys.push_back(static_cast<int>(MOD));

    std::vector<Event> events;
    events.reserve(2 * t);

    for (const auto& r : rects) {
        ys.push_back(r.y1);
        ys.push_back(r.y2 + 1);
        events.push_back({r.x1, r.y1, r.y2, +1});
        if (r.x2 + 1 < MOD) {
            events.push_back({r.x2 + 1, r.y1, r.y2, -1});
        }
    }

    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

    std::sort(events.begin(), events.end());

    SegTree seg(ys);

    u64 cnt[12] = {};

    int cur_x = 0;
    std::size_t ei = 0;

    while (ei < events.size()) {
        const int nx = events[ei].x;
        const i64 width = static_cast<i64>(nx) - static_cast<i64>(cur_x);
        if (width > 0) {
            const auto& a = seg.all();
            for (int r = 0; r < 12; ++r) {
                cnt[r] += static_cast<u64>(a[r]) * static_cast<u64>(width);
            }
        }

        while (ei < events.size() && events[ei].x == nx) {
            const auto& e = events[ei];
            const int l = static_cast<int>(std::lower_bound(ys.begin(), ys.end(), e.y1) - ys.begin());
            const int rr = static_cast<int>(std::lower_bound(ys.begin(), ys.end(), e.y2 + 1) - ys.begin()) - 1;
            seg.update(l, rr, e.delta);
            ++ei;
        }

        cur_x = nx;
    }

    if (cur_x < MOD) {
        const i64 width = static_cast<i64>(MOD) - static_cast<i64>(cur_x);
        const auto& a = seg.all();
        for (int r = 0; r < 12; ++r) {
            cnt[r] += static_cast<u64>(a[r]) * static_cast<u64>(width);
        }
    }

    u64 ans = 0;
    ans += cnt[0] * 12ULL;
    for (int r = 1; r <= 11; ++r) {
        ans += cnt[r] * static_cast<u64>(r);
    }

    return ans;
}

int main() {
    assert(solve(0) == 30'621'295'449'583'788ULL);
    assert(solve(1) == 30'613'048'345'941'659ULL);
    assert(solve(10) == 21'808'930'308'198'471ULL);
    assert(solve(100) == 16'190'667'393'984'172ULL);

    std::cout << solve(100'000) << '\n';
    return 0;
}
