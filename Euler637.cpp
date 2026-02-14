#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;

struct BitMin {
    int n;
    std::vector<int> bit;
    std::vector<int> arr;
    static constexpr int INF = 1'000'000'000;

    explicit BitMin(const std::vector<int>& initial)
        : n(static_cast<int>(initial.size()) - 1), bit(static_cast<std::size_t>(n + 1), INF), arr(initial) {
        for (int i = 1; i <= n; ++i) update_bit_min(rev_index(i), arr[static_cast<std::size_t>(i)]);
    }

    void update(int i, int new_val) {
        if (new_val >= arr[static_cast<std::size_t>(i)]) return;
        arr[static_cast<std::size_t>(i)] = new_val;
        update_bit_min(rev_index(i), new_val);
    }

    int suffix_min(int i) const {
        if (i == 0) return 0;
        int ans = arr[static_cast<std::size_t>(i)];
        for (int x = rev_index(i); x > 0; x -= x & -x) ans = std::min(ans, bit[static_cast<std::size_t>(x)]);
        return ans;
    }

    int rev_index(int idx) const { return n - idx + 1; }

    void update_bit_min(int idx, int val) {
        for (int x = idx; x <= n; x += x & -x) {
            int& cur = bit[static_cast<std::size_t>(x)];
            if (val < cur) cur = val;
        }
    }
};

static int digit_sum(int n, int b) {
    int s = 0;
    while (n > 0) {
        s += n % b;
        n /= b;
    }
    return s;
}

struct Node {
    int krem;
    int digit_sum_rem;
    int exp;
    int curr_part;
    int part_sum;
};

static std::vector<int> f_steps(int n, int b) {
    std::vector<int> steps(static_cast<std::size_t>(n + 1), 0);
    for (int i = b; i <= n; ++i) steps[static_cast<std::size_t>(i)] = n;

    BitMin bit_min(steps);

    std::vector<int> b_pows;
    for (i64 p = 1; p <= n; p *= b) b_pows.push_back(static_cast<int>(p));

    std::vector<Node> st;
    st.reserve(1 << 14);

    for (int k = b; k <= n; ++k) {
        const int dsum = digit_sum(k, b);
        int current_min = steps[static_cast<std::size_t>(dsum)] + 1;

        st.clear();
        st.push_back({k, dsum, 0, 0, 0});

        while (!st.empty()) {
            const Node cur = st.back();
            st.pop_back();

            const int smallest = cur.digit_sum_rem + cur.part_sum + cur.curr_part;
            const int cand = steps[static_cast<std::size_t>(smallest)] + 1;
            if (cand < current_min) current_min = cand;

            if (bit_min.suffix_min(smallest) + 1 >= current_min) continue;
            if (cur.krem == 0) continue;

            const int digit = cur.krem % b;
            const int next_k = cur.krem / b;
            const int next_sum = cur.digit_sum_rem - digit;

            if (cur.exp != 0) {
                const int next_curr = cur.curr_part + digit * b_pows[static_cast<std::size_t>(cur.exp)];
                st.push_back({next_k, next_sum, cur.exp + 1, next_curr, cur.part_sum});
            }
            st.push_back({next_k, next_sum, 1, digit, cur.part_sum + cur.curr_part});
        }

        steps[static_cast<std::size_t>(k)] = current_min;
        bit_min.update(k, current_min);
    }

    return steps;
}

static i64 g(int n, int b1, int b2) {
    std::vector<int> s1 = f_steps(n, b1);
    std::vector<int> s2 = f_steps(n, b2);
    i64 total = 0;
    for (int k = 1; k <= n; ++k) {
        if (s1[static_cast<std::size_t>(k)] == s2[static_cast<std::size_t>(k)]) total += k;
    }
    return total;
}

int main() {
    assert(g(100, 10, 3) == 3302);
    std::cout << g(10'000'000, 10, 3) << '\n';
    return 0;
}
