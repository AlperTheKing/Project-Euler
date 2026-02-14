#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

constexpr int MOD = 998244353;

inline int add_mod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

inline int sub_mod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

inline int mul_mod(long long a, long long b) {
    return static_cast<int>((static_cast<__int128>(a) * b) % MOD);
}

enum class EdgeType : std::uint8_t { G = 0, O = 1, C = 2 };

struct Edge {
    EdgeType t{EdgeType::G};
    int v{0};
};

inline Edge G() { return Edge{EdgeType::G, 0}; }
inline Edge O() { return Edge{EdgeType::O, 0}; }
inline Edge C(int n) { return Edge{EdgeType::C, n}; }

inline bool is_g(const Edge& e) { return e.t == EdgeType::G; }
inline bool is_o(const Edge& e) { return e.t == EdgeType::O; }
inline bool is_c(const Edge& e) { return e.t == EdgeType::C; }

Edge insert_edge(int n, const Edge& e) {
    if (is_g(e)) return C(n);
    if (is_o(e)) return (n & 1) ? C(n) : C(0);
    assert((e.v & n) > 0);
    return C(n);
}

Edge up0(const Edge& e) {
    if (is_c(e)) return C(e.v / 2);
    if (is_o(e)) return C(0);
    return G();
}

Edge up1(const Edge& e) {
    if (is_c(e)) {
        if ((e.v & 1) == 0) return C(e.v / 2);
        return G();
    }
    return G();
}

struct Key {
    int l;
    int n;
    int left_code;
    int right_code;

    bool operator==(const Key& other) const noexcept {
        return l == other.l && n == other.n && left_code == other.left_code &&
               right_code == other.right_code;
    }
};

struct KeyHash {
    std::size_t operator()(const Key& k) const noexcept {
        std::size_t h = std::hash<int>{}(k.l);
        auto mix = [&](int x) {
            std::size_t v = std::hash<int>{}(x);
            h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
        };
        mix(k.n);
        mix(k.left_code);
        mix(k.right_code);
        return h;
    }
};

int edge_code(const Edge& e) {
    if (is_g(e)) return -1;
    if (is_o(e)) return -2;
    return e.v;
}

class Solver774 {
  public:
    Solver774(int max_l, int a, int b) : a_(a), b_(b) {
        fib_pos_.assign(max_l + 2, 0);
        fib_pos_[0] = 0;
        fib_pos_[1] = 1;
        for (int i = 2; i < static_cast<int>(fib_pos_.size()); ++i) {
            fib_pos_[i] = add_mod(fib_pos_[i - 1], fib_pos_[i - 2]);
        }
        memo_.reserve(1 << 20U);
    }

    int solve() { return f(a_, b_, G(), G()); }

  private:
    int a_;
    int b_;
    std::vector<int> fib_pos_;
    std::unordered_map<Key, int, KeyHash> memo_;

    int fib(int x) const {
        if (x >= 0) return fib_pos_[static_cast<std::size_t>(x)];
        int n = -x;
        int v = fib_pos_[static_cast<std::size_t>(n)];
        if ((n & 1) == 0 && v != 0) v = MOD - v;
        return v;
    }

    int f(int l, int n, const Edge& left, const Edge& right) {
        if ((is_c(left) && left.v == 0) || (is_c(right) && right.v == 0)) return 0;

        const Key key{l, n, edge_code(left), edge_code(right)};
        auto it = memo_.find(key);
        if (it != memo_.end()) return it->second;

        int ret = 0;

        if (n == 0) {
            ret = (l <= 1 && is_g(left) && is_g(right)) ? 1 : 0;
            memo_.emplace(key, ret);
            return ret;
        }

        if (n == 1) {
            if (is_g(left) && is_g(right)) {
                ret = (l > 1) ? 1 : 2;
            } else if ((is_o(left) && is_o(right)) || (is_o(left) && is_g(right)) ||
                       (is_g(left) && is_o(right))) {
                ret = 1;
            } else if (is_c(left) && is_c(right)) {
                ret = ((left.v & 1) && (right.v & 1)) ? 1 : 0;
            } else if (is_c(left)) {
                ret = (left.v & 1) ? 1 : 0;
            } else if (is_c(right)) {
                ret = (right.v & 1) ? 1 : 0;
            } else {
                ret = 0;
            }
            memo_.emplace(key, ret);
            return ret;
        }

        const int m = (n - 1) / 2;

        if (l == 1) {
            if (is_g(left) && is_g(right)) {
                ret = (n + 1) % MOD;
            } else {
                ret = add_mod(f(1, n / 2, up0(left), up0(right)),
                              f(1, m, up1(left), up1(right)));
            }
            memo_.emplace(key, ret);
            return ret;
        }

        int s = 0;
        s = add_mod(s, mul_mod(f(l, m, up0(left), up0(right)), fib(l)));
        s = add_mod(s, mul_mod(f(l, m, up0(left), up1(right)), fib(l - 1)));
        s = add_mod(s, mul_mod(f(l, m, up1(left), up0(right)), fib(l - 1)));
        s = add_mod(s, mul_mod(f(l, m, up1(left), up1(right)), fib(l - 2)));

        const Edge up1o = up1(O());
        for (int x = 1; x <= l - 1; ++x) {
            int right_part = f(l - x, n, O(), right);
            if (x > 1) {
                int left_part = f(x, m, up0(left), up1o);
                s = add_mod(s, mul_mod(mul_mod(left_part, right_part), fib(x - 1)));
            }
            int left_part2 = f(x, m, up1(left), up1o);
            s = add_mod(s, mul_mod(mul_mod(left_part2, right_part), fib(x - 2)));
        }

        if ((n & 1) == 0) {
            const Edge cn = C(n);
            const Edge up0cn = up0(cn);

            for (int x = 2; x <= l - 1; ++x) {
                int right_part = f(l - x, n, cn, right);

                int a = f(x - 1, m, up0(left), up0cn);
                s = add_mod(s, mul_mod(mul_mod(a, right_part), fib(x)));

                int b = f(x - 1, m, up1(left), up0cn);
                s = add_mod(s, mul_mod(mul_mod(b, right_part), fib(x - 1)));
            }

            s = add_mod(s, f(l - 1, n, insert_edge(n, left), right));

            const Edge ins_right = insert_edge(n, right);
            s = add_mod(s, mul_mod(f(l - 1, m, up0(left), up0(ins_right)), fib(l)));
            s = add_mod(s, mul_mod(f(l - 1, m, up1(left), up0(ins_right)), fib(l - 1)));
        }

        ret = s;
        memo_.emplace(key, ret);
        return ret;
    }
};

int compute_c(int n, int b) {
    Solver774 solver(n, n, b);
    return solver.solve();
}

}  // namespace

int main() {
    struct Check {
        int n;
        std::uint64_t b;
        int expected;
    };

    const Check checks[] = {
        {3, 4, 18},
        {10, 6, 2496120},
        {100, 200, 268159379},
    };

    for (const auto& chk : checks) {
        int got = compute_c(chk.n, static_cast<int>(chk.b));
        if (got != chk.expected) {
            std::cerr << "Validation failure: c(" << chk.n << ", " << chk.b
                      << ") = " << got << ", expected " << chk.expected << '\n';
            return 1;
        }
    }

    std::cout << compute_c(123, 123456789) << '\n';
    return 0;
}
