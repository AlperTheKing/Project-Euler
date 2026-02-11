#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

struct State {
    int next[2];
    int link;
    int len;
    int occ;
};

std::vector<int> build_L(int n) {
    std::vector<State> st(static_cast<std::size_t>(2 * n));

    int sz = 1;
    int last = 0;

    st[0].next[0] = st[0].next[1] = -1;
    st[0].link = -1;
    st[0].len = 0;
    st[0].occ = 0;

    const long double invphi = (std::sqrt(5.0L) - 1.0L) / 2.0L;
    std::uint64_t beat_prev = 0ULL;

    for (int i = 0; i < n; ++i) {
        const int a = __builtin_popcount(static_cast<unsigned int>(i)) & 1;

        const std::uint64_t beat_cur =
            static_cast<std::uint64_t>(std::floor((static_cast<long double>(i) + 1.0L) * invphi));
        const int b = static_cast<int>(beat_cur - beat_prev);
        beat_prev = beat_cur;

        const int c = a ^ b;

        const int cur = sz++;
        st[cur].next[0] = st[cur].next[1] = -1;
        st[cur].len = st[last].len + 1;
        st[cur].occ = 1;

        int p = last;
        while (p != -1 && st[p].next[c] == -1) {
            st[p].next[c] = cur;
            p = st[p].link;
        }

        if (p == -1) {
            st[cur].link = 0;
        } else {
            const int q = st[p].next[c];
            if (st[p].len + 1 == st[q].len) {
                st[cur].link = q;
            } else {
                const int clone = sz++;
                st[clone] = st[q];
                st[clone].len = st[p].len + 1;
                st[clone].occ = 0;

                while (p != -1 && st[p].next[c] == q) {
                    st[p].next[c] = clone;
                    p = st[p].link;
                }

                st[q].link = clone;
                st[cur].link = clone;
            }
        }

        last = cur;
    }

    st.resize(static_cast<std::size_t>(sz));

    std::vector<int> cnt_len(static_cast<std::size_t>(n + 1), 0);
    for (int i = 0; i < sz; ++i) {
        ++cnt_len[static_cast<std::size_t>(st[i].len)];
    }
    for (int i = 1; i <= n; ++i) {
        cnt_len[static_cast<std::size_t>(i)] += cnt_len[static_cast<std::size_t>(i - 1)];
    }

    std::vector<int> order(static_cast<std::size_t>(sz), 0);
    for (int i = sz - 1; i >= 0; --i) {
        const int l = st[i].len;
        order[static_cast<std::size_t>(--cnt_len[static_cast<std::size_t>(l)])] = i;
    }

    for (int i = sz - 1; i > 0; --i) {
        const int v = order[static_cast<std::size_t>(i)];
        const int p = st[v].link;
        if (p >= 0) {
            st[p].occ += st[v].occ;
        }
    }

    std::vector<int> best(static_cast<std::size_t>(n + 1), 0);
    for (int i = 1; i < sz; ++i) {
        const int occ = st[i].occ;
        if (occ <= n && st[i].len > best[static_cast<std::size_t>(occ)]) {
            best[static_cast<std::size_t>(occ)] = st[i].len;
        }
    }

    for (int k = n - 1; k >= 1; --k) {
        if (best[static_cast<std::size_t>(k)] < best[static_cast<std::size_t>(k + 1)]) {
            best[static_cast<std::size_t>(k)] = best[static_cast<std::size_t>(k + 1)];
        }
    }

    return best;
}

std::uint64_t sum_non_zero_L(int n) {
    const std::vector<int> L = build_L(n);
    std::uint64_t out = 0ULL;
    for (int k = 1; k <= n; ++k) {
        const int v = L[static_cast<std::size_t>(k)];
        if (v == 0) {
            break;
        }
        out += static_cast<std::uint64_t>(v);
    }
    return out;
}

}  // namespace

int main() {
    {
        const std::vector<int> L10 = build_L(10);
        assert(L10[2] == 5);
        assert(L10[3] == 2);
    }

    {
        const std::vector<int> L100 = build_L(100);
        assert(L100[2] == 14);
        assert(L100[4] == 6);
    }

    {
        const std::vector<int> L1000 = build_L(1000);
        assert(L1000[2] == 86);
        assert(L1000[3] == 45);
        assert(L1000[5] == 31);
    }

    assert(sum_non_zero_L(1000) == 2460ULL);

    std::cout << sum_non_zero_L(5'000'000) << "\n";
    return 0;
}
