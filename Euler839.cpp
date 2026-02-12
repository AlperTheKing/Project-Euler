#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Block {
    u64 sum;
    u32 len;
};

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int d = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u128 compute_B(u32 n) {
    static constexpr u64 kMod = 50'515'093ULL;

    std::vector<Block> st;
    st.reserve(1'000'000);

    u64 s = 290'797ULL;
    u128 pref_init = 0;
    u128 sum_pref_init = 0;

    for (u32 i = 0; i < n; ++i) {
        const u64 v = s;
        pref_init += static_cast<u128>(v);
        sum_pref_init += pref_init;

        s = static_cast<u64>((static_cast<u128>(s) * s) % kMod);

        st.push_back({v, 1});
        while (st.size() >= 2) {
            Block& right = st.back();
            Block& left = st[st.size() - 2];
            const u64 left_max = (left.sum + left.len - 1ULL) / left.len;
            const u64 right_min = right.sum / right.len;
            if (left_max <= right_min) break;
            left.sum += right.sum;
            left.len += right.len;
            st.pop_back();
        }
    }

    u128 pref_final = 0;
    u128 sum_pref_final = 0;

    for (const Block& b : st) {
        const u64 len = b.len;
        const u64 q = b.sum / len;
        const u64 r = b.sum % len;

        const u64 cnt_q = len - r;
        if (cnt_q > 0) {
            sum_pref_final += static_cast<u128>(cnt_q) * pref_final;
            sum_pref_final += static_cast<u128>(q) * cnt_q * (cnt_q + 1ULL) / 2ULL;
            pref_final += static_cast<u128>(cnt_q) * q;
        }

        const u64 cnt_q1 = r;
        if (cnt_q1 > 0) {
            const u64 qp1 = q + 1ULL;
            sum_pref_final += static_cast<u128>(cnt_q1) * pref_final;
            sum_pref_final += static_cast<u128>(qp1) * cnt_q1 * (cnt_q1 + 1ULL) / 2ULL;
            pref_final += static_cast<u128>(cnt_q1) * qp1;
        }
    }

    assert(pref_init == pref_final);
    return sum_pref_init - sum_pref_final;
}

int main() {
    assert(compute_B(5) == 0);
    assert(compute_B(6) == 14'263'289ULL);
    assert(compute_B(100) == 3'284'417'556ULL);

    std::cout << to_string_u128(compute_B(10'000'000)) << '\n';
    return 0;
}
