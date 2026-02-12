#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

struct Node {
    int child[2];
    std::uint32_t cnt_term;
    Node() : child{-1, -1}, cnt_term(0) {}
};

class PrimeGuessing {
public:
    explicit PrimeGuessing(int limit) : limit_(limit) {
        nodes_.reserve(32'000'000);
        nodes_.push_back(Node());
        build();
    }

    long double expectation() {
        return dfs(0);
    }

private:
    static constexpr std::uint32_t TERM_BIT = 0x80000000u;
    static constexpr std::uint32_t CNT_MASK = 0x7fffffffu;

    int limit_;
    std::vector<Node> nodes_;

    static inline std::uint32_t count_of(const Node& n) {
        return n.cnt_term & CNT_MASK;
    }

    static inline std::uint32_t term_of(const Node& n) {
        return n.cnt_term >> 31;
    }

    inline void inc_count(int idx) {
        Node& n = nodes_[idx];
        std::uint32_t c = (n.cnt_term & CNT_MASK) + 1u;
        n.cnt_term = (n.cnt_term & TERM_BIT) | c;
    }

    inline void set_term(int idx) {
        nodes_[idx].cnt_term |= TERM_BIT;
    }

    void insert_prime(std::uint32_t p) {
        int u = 0;
        inc_count(u);

        int len = 32 - __builtin_clz(p);
        for (int d = 0; d < len; ++d) {
            int b = (p >> d) & 1u;
            int v = nodes_[u].child[b];
            if (v == -1) {
                v = static_cast<int>(nodes_.size());
                nodes_[u].child[b] = v;
                nodes_.push_back(Node());
            }
            u = v;
            inc_count(u);
        }
        set_term(u);
    }

    void build() {
        if (limit_ >= 2) insert_prime(2);
        if (limit_ < 3) return;

        int m = (limit_ >> 1) + 1;
        std::vector<std::uint8_t> is_prime(m, 1);
        is_prime[0] = 0;

        int r = static_cast<int>(std::sqrt(static_cast<long double>(limit_)));
        for (int p = 3; p <= r; p += 2) {
            if (!is_prime[p >> 1]) continue;
            int step = p << 1;
            for (int x = p * p; x <= limit_; x += step) {
                is_prime[x >> 1] = 0;
            }
        }

        for (int p = 3; p <= limit_; p += 2) {
            if (is_prime[p >> 1]) insert_prime(static_cast<std::uint32_t>(p));
        }
    }

    long double dfs(int u) {
        const Node& cur = nodes_[u];
        std::uint32_t total = count_of(cur);
        std::uint32_t term = term_of(cur);
        std::uint32_t cont = total - term;
        if (cont == 0) return 0.0L;

        std::uint32_t c0 = 0, c1 = 0;
        std::uint32_t n0 = 0, n1 = 0;
        long double v0 = 0.0L, v1 = 0.0L;

        int ch0 = cur.child[0];
        if (ch0 != -1) {
            c0 = count_of(nodes_[ch0]);
            n0 = c0 - term_of(nodes_[ch0]);
            if (n0 > 0) v0 = dfs(ch0);
        }

        int ch1 = cur.child[1];
        if (ch1 != -1) {
            c1 = count_of(nodes_[ch1]);
            n1 = c1 - term_of(nodes_[ch1]);
            if (n1 > 0) v1 = dfs(ch1);
        }

        long double denom = static_cast<long double>(cont);
        long double score = static_cast<long double>(std::max(c0, c1)) / denom;
        score += static_cast<long double>(n0) / denom * v0;
        score += static_cast<long double>(n1) / denom * v1;
        return score;
    }
};

int main() {
    {
        PrimeGuessing g(10);
        assert(std::fabsl(g.expectation() - 2.0L) < 1e-15L);
    }
    {
        PrimeGuessing g(30);
        assert(std::fabsl(g.expectation() - 2.9L) < 1e-15L);
    }

    PrimeGuessing g(100'000'000);
    std::cout << std::fixed << std::setprecision(8) << g.expectation() << '\n';
    return 0;
}
