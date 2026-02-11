#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr int kSeqMod = 10'000'019;
constexpr int kAnswerMod = 1'000'000'007;

std::vector<int> generate_sequence(const int n) {
    std::vector<int> a(static_cast<std::size_t>(n));
    i64 cur = 1;
    for (int i = 0; i < n; ++i) {
        cur = (cur * 153LL) % kSeqMod;
        a[static_cast<std::size_t>(i)] = static_cast<int>(cur);
    }
    return a;
}

u64 brute_exact_sum4(const std::vector<int>& a) {
    const int n = static_cast<int>(a.size());
    u64 total = 0ULL;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (a[static_cast<std::size_t>(i)] >= a[static_cast<std::size_t>(j)]) {
                continue;
            }
            for (int k = j + 1; k < n; ++k) {
                if (a[static_cast<std::size_t>(j)] >= a[static_cast<std::size_t>(k)]) {
                    continue;
                }
                for (int l = k + 1; l < n; ++l) {
                    if (a[static_cast<std::size_t>(k)] >= a[static_cast<std::size_t>(l)]) {
                        continue;
                    }
                    total += static_cast<u64>(a[static_cast<std::size_t>(i)]) +
                             static_cast<u64>(a[static_cast<std::size_t>(j)]) +
                             static_cast<u64>(a[static_cast<std::size_t>(k)]) +
                             static_cast<u64>(a[static_cast<std::size_t>(l)]);
                }
            }
        }
    }
    return total;
}

struct FenwickPair {
    int n = 0;
    std::vector<int> cnt;
    std::vector<int> sum;

    explicit FenwickPair(const int n_) : n(n_), cnt(static_cast<std::size_t>(n_ + 1), 0), sum(static_cast<std::size_t>(n_ + 1), 0) {}

    static int add_mod(const int a, const int b) {
        int s = a + b;
        if (s >= kAnswerMod) {
            s -= kAnswerMod;
        }
        return s;
    }

    std::pair<int, int> query(int idx) const {
        int c = 0;
        int s = 0;
        while (idx > 0) {
            c += cnt[static_cast<std::size_t>(idx)];
            if (c >= kAnswerMod) {
                c -= kAnswerMod;
            }
            s += sum[static_cast<std::size_t>(idx)];
            if (s >= kAnswerMod) {
                s -= kAnswerMod;
            }
            idx -= idx & -idx;
        }
        return {c, s};
    }

    void update(int idx, const int add_cnt, const int add_sum) {
        while (idx <= n) {
            cnt[static_cast<std::size_t>(idx)] = add_mod(cnt[static_cast<std::size_t>(idx)], add_cnt);
            sum[static_cast<std::size_t>(idx)] = add_mod(sum[static_cast<std::size_t>(idx)], add_sum);
            idx += idx & -idx;
        }
    }
};

int solve_mod(const int n) {
    const std::vector<int> a = generate_sequence(n);

    std::vector<int> sorted = a;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

    std::vector<int> rank(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        rank[static_cast<std::size_t>(i)] =
            static_cast<int>(std::lower_bound(sorted.begin(), sorted.end(), a[static_cast<std::size_t>(i)]) - sorted.begin()) +
            1;
    }

    const int m = static_cast<int>(sorted.size());
    std::array<FenwickPair, 3> fw = {FenwickPair(m), FenwickPair(m), FenwickPair(m)};

    int answer = 0;

    for (int i = 0; i < n; ++i) {
        const int r = rank[static_cast<std::size_t>(i)];
        const int x = a[static_cast<std::size_t>(i)];

        const int c1 = 1;
        const int s1 = x;

        const auto [c1_pref, s1_pref] = fw[0].query(r - 1);
        const int c2 = c1_pref;
        const int s2 = static_cast<int>((static_cast<i64>(s1_pref) + static_cast<i64>(x) * c2) % kAnswerMod);

        const auto [c2_pref, s2_pref] = fw[1].query(r - 1);
        const int c3 = c2_pref;
        const int s3 = static_cast<int>((static_cast<i64>(s2_pref) + static_cast<i64>(x) * c3) % kAnswerMod);

        const auto [c3_pref, s3_pref] = fw[2].query(r - 1);
        const int c4 = c3_pref;
        const int s4 = static_cast<int>((static_cast<i64>(s3_pref) + static_cast<i64>(x) * c4) % kAnswerMod);

        answer += s4;
        if (answer >= kAnswerMod) {
            answer -= kAnswerMod;
        }

        fw[0].update(r, c1, s1);
        fw[1].update(r, c2, s2);
        fw[2].update(r, c3, s3);
    }

    return answer;
}

}  // namespace

int main() {
    {
        const std::vector<int> a6 = generate_sequence(6);
        assert(brute_exact_sum4(a6) == 94'513'710ULL);
    }
    {
        const std::vector<int> a100 = generate_sequence(100);
        assert(brute_exact_sum4(a100) == 4'465'488'724'217ULL);
    }

    std::cout << solve_mod(1'000'000) << '\n';
    return 0;
}
