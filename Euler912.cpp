#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1'000'000'007LL;

struct Stats {
    u64 total = 0;
    u64 odd = 0;
    i64 odd_mod = 0;
    i64 s1_mod = 0;
    i64 s2_mod = 0;
};

i64 add_mod(i64 a, i64 b) {
    i64 v = a + b;
    if (v >= kMod) {
        v -= kMod;
    }
    return v;
}

i64 mul_mod(i64 a, i64 b) {
    return static_cast<i64>((__int128)a * b % kMod);
}

class Solver {
public:
    Solver() {
        for (int c = 0; c < 3; ++c) {
            memo_[c].resize(1);
            seen_[c].resize(1, false);
        }
    }

    i64 F(u64 n) {
        u64 remaining = n;
        u64 prefix = 0;
        i64 ans = 0;
        int length = 1;

        while (remaining > 0) {
            const Stats block = full(1, length - 1);
            const u64 take = (remaining < block.total) ? remaining : block.total;
            const Stats chosen = (take == block.total) ? block : prefix_stats(1, length - 1, take);

            const i64 p = static_cast<i64>(prefix % kMod);
            const i64 p2 = mul_mod(p, p);

            i64 contrib = 0;
            contrib = add_mod(contrib, mul_mod(chosen.odd_mod, p2));
            contrib = add_mod(contrib, mul_mod((2 * p) % kMod, chosen.s1_mod));
            contrib = add_mod(contrib, chosen.s2_mod);

            ans = add_mod(ans, contrib);
            prefix += take;
            remaining -= take;
            ++length;
        }

        return ans;
    }

private:
    std::array<std::vector<Stats>, 3> memo_;
    std::array<std::vector<char>, 3> seen_;

    Stats empty() const {
        return Stats{};
    }

    Stats merge(const Stats& left, const Stats& right) const {
        if (left.total == 0) {
            return right;
        }
        if (right.total == 0) {
            return left;
        }

        Stats out;
        out.total = left.total + right.total;
        out.odd = left.odd + right.odd;
        out.odd_mod = add_mod(left.odd_mod, right.odd_mod);

        const i64 shift = static_cast<i64>(left.total % kMod);
        const i64 shift2 = mul_mod(shift, shift);

        out.s1_mod = left.s1_mod;
        out.s1_mod = add_mod(out.s1_mod, right.s1_mod);
        out.s1_mod = add_mod(out.s1_mod, mul_mod(shift, right.odd_mod));

        out.s2_mod = left.s2_mod;
        out.s2_mod = add_mod(out.s2_mod, right.s2_mod);
        out.s2_mod = add_mod(out.s2_mod, mul_mod((2 * shift) % kMod, right.s1_mod));
        out.s2_mod = add_mod(out.s2_mod, mul_mod(shift2, right.odd_mod));

        return out;
    }

    Stats base(int c) const {
        Stats s;
        s.total = 1;
        s.odd = (c > 0) ? 1 : 0;
        s.odd_mod = static_cast<i64>(s.odd);
        s.s1_mod = static_cast<i64>(s.odd);
        s.s2_mod = static_cast<i64>(s.odd);
        return s;
    }

    Stats full(int c, int r) {
        if (r >= static_cast<int>(memo_[c].size())) {
            memo_[c].resize(static_cast<std::size_t>(r + 1));
            seen_[c].resize(static_cast<std::size_t>(r + 1), false);
        }

        if (seen_[c][r]) {
            return memo_[c][r];
        }

        Stats out;
        if (r == 0) {
            out = base(c);
        } else {
            const Stats left = full(0, r - 1);
            const Stats right = (c < 2) ? full(c + 1, r - 1) : empty();
            out = merge(left, right);
        }

        seen_[c][r] = true;
        memo_[c][r] = out;
        return out;
    }

    Stats prefix_stats(int c, int r, u64 k) {
        if (k == 0) {
            return empty();
        }

        const Stats f = full(c, r);
        if (k >= f.total) {
            return f;
        }

        if (r == 0) {
            return base(c);
        }

        const Stats left = full(0, r - 1);
        if (k <= left.total) {
            return prefix_stats(0, r - 1, k);
        }

        assert(c < 2);
        const Stats right_part = prefix_stats(c + 1, r - 1, k - left.total);
        return merge(left, right_part);
    }
};

bool no_three_ones(u64 x) {
    int run = 0;
    while (x > 0) {
        if (x & 1ULL) {
            ++run;
            if (run == 3) {
                return false;
            }
        } else {
            run = 0;
        }
        x >>= 1;
    }
    return true;
}

i64 brute_F(int n) {
    i64 ans = 0;
    int idx = 0;
    for (u64 x = 1; idx < n; ++x) {
        if (!no_three_ones(x)) {
            continue;
        }
        ++idx;
        if (x & 1ULL) {
            ans += 1LL * idx * idx;
            ans %= kMod;
        }
    }
    return ans;
}

void validate() {
    Solver solver;
    assert(solver.F(10) == 199);

    for (int n = 1; n <= 200; ++n) {
        assert(solver.F(static_cast<u64>(n)) == brute_F(n));
    }
}

}  // namespace

int main() {
    validate();
    Solver solver;
    std::cout << solver.F(10'000'000'000'000'000ULL) << '\n';
    return 0;
}
