#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    i64 n = 1000000000000000000LL;
    i64 k = 1000000;
    i64 mod = 1000000000;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) return false;
    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') return false;
        parsed = parsed * 10 + (c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_i64_after_prefix(arg, "--N=", options.n)) continue;
        if (parse_i64_after_prefix(arg, "--K=", options.k)) continue;
        if (parse_i64_after_prefix(arg, "--mod=", options.mod)) continue;
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n > 0 && options.k >= 0 && options.mod > 0;
}

inline i64 add_mod(i64 a, i64 b, i64 mod) {
    const i64 threshold = mod - b;
    return a >= threshold ? a - threshold : a + b;
}

inline i64 sub_mod(i64 a, i64 b, i64 mod) {
    return a >= b ? a - b : mod - (b - a);
}

inline i64 mul_mod(i64 a, i64 b, i64 mod) {
    return static_cast<i64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

inline i64 add_mod_n(i64 a, i64 b, i64 n) {
    const i64 threshold = n - b;
    return a >= threshold ? a - threshold : a + b;
}

struct Node {
    int parent = 0;
    int child[2]{0, 0};
    bool rev = false;
    i64 s = 0;
    i64 e = 0;
    i64 cnt = 0;
};

struct SplayTree {
    i64 n = 0;
    i64 mod = 1;
    std::vector<Node> nodes;

    explicit SplayTree(i64 n_, i64 mod_) : n(n_), mod(mod_) {
        nodes.resize(3000500);
        init();
    }

    void init() {
        nodes[0].child[0] = 1;
        nodes[0].cnt = n + 2;

        nodes[1].parent = 0;
        nodes[1].s = 0;
        nodes[1].e = n - 1;
        nodes[1].child[0] = 2;
        nodes[1].child[1] = 3;
        nodes[1].cnt = n + 2;

        nodes[2].parent = 1;
        nodes[2].s = -1;
        nodes[2].e = -1;
        nodes[2].cnt = 1;

        nodes[3].parent = 1;
        nodes[3].s = n;
        nodes[3].e = n;
        nodes[3].cnt = 1;

        top = 10;
    }

    int top = 0;

    static i64 seg_len(const Node& nd) {
        return std::llabs(nd.s - nd.e) + 1;
    }

    void push_down(int x) {
        Node& nd = nodes[x];
        if (!nd.rev) return;
        nd.rev = false;
        for (int dir = 0; dir < 2; ++dir) {
            if (nd.child[dir]) {
                nodes[nd.child[dir]].rev = !nodes[nd.child[dir]].rev;
                std::swap(nodes[nd.child[dir]].s, nodes[nd.child[dir]].e);
            }
        }
        std::swap(nd.child[0], nd.child[1]);
    }

    void update(int x) {
        Node& nd = nodes[x];
        i64 cnt = seg_len(nd);
        if (nd.child[0]) cnt += nodes[nd.child[0]].cnt;
        if (nd.child[1]) cnt += nodes[nd.child[1]].cnt;
        nd.cnt = cnt;
    }

    void rotate(int x, int dir) {
        int y = nodes[x].parent;
        push_down(y);
        push_down(x);

        nodes[y].child[dir ^ 1] = nodes[x].child[dir];
        if (nodes[x].child[dir]) nodes[nodes[x].child[dir]].parent = y;

        nodes[x].parent = nodes[y].parent;
        if (nodes[nodes[y].parent].child[0] == y) {
            nodes[nodes[y].parent].child[0] = x;
        } else {
            nodes[nodes[y].parent].child[1] = x;
        }

        nodes[x].child[dir] = y;
        nodes[y].parent = x;
        update(y);
    }

    void splay(int x, int target_parent) {
        for (push_down(x); nodes[x].parent != target_parent;) {
            int y = nodes[x].parent;
            int z = nodes[y].parent;
            if (z == target_parent) {
                if (nodes[y].child[0] == x) rotate(x, 1);
                else rotate(x, 0);
            } else {
                if (nodes[z].child[0] == y) {
                    if (nodes[y].child[0] == x) {
                        rotate(y, 1);
                        rotate(x, 1);
                    } else {
                        rotate(x, 0);
                        rotate(x, 1);
                    }
                } else {
                    if (nodes[y].child[1] == x) {
                        rotate(y, 0);
                        rotate(x, 0);
                    } else {
                        rotate(x, 1);
                        rotate(x, 0);
                    }
                }
            }
        }
        update(x);
    }

    int find(i64& key) {
        int curr = nodes[0].child[0];
        while (curr) {
            push_down(curr);
            i64 left_size = nodes[curr].child[0] ? nodes[nodes[curr].child[0]].cnt : 0;
            if (key <= left_size) {
                curr = nodes[curr].child[0];
                continue;
            }
            key -= left_size;
            i64 mid = seg_len(nodes[curr]);
            if (key <= mid) return curr;
            key -= mid;
            curr = nodes[curr].child[1];
        }
        return curr;
    }

    i64 cal(i64 s, i64 e, i64 idx) const {
        const i64 delta = std::llabs(s - e);
        const i64 mod2 = mod * 2;
        const i64 mod6 = mod * 6;
        const i64 is1 = static_cast<i64>((static_cast<u128>(delta) % mod2) * (delta + 1) % mod2 / 2 % mod);
        const i64 is2 = static_cast<i64>((static_cast<u128>(delta) % mod6) * (delta + 1) % mod6 *
                                         ((2 * delta + 1) % mod6) % mod6 / 6 % mod);

        if (s <= e) {
            i64 t = mul_mod(idx % mod, s % mod, mod);
            t = mul_mod(t, (delta + 1) % mod, mod);
            i64 term = mul_mod((idx + s) % mod, is1, mod);
            t = add_mod(t, term, mod);
            t = add_mod(t, is2 % mod, mod);
            return t;
        }
        i64 t = mul_mod(idx % mod, s % mod, mod);
        t = mul_mod(t, (delta + 1) % mod, mod);
        i64 term = mul_mod(((s - idx) % mod + mod) % mod, is1, mod);
        t = add_mod(t, term, mod);
        t = sub_mod(t, is2 % mod, mod);
        return t;
    }

    i64 idx = 0;
    i64 solve_impl(int r) {
        i64 ret = 0;
        push_down(r);
        if (nodes[r].child[0]) ret = add_mod(ret, solve_impl(nodes[r].child[0]), mod);

        const i64 len = seg_len(nodes[r]);
        if (idx >= 0 && idx + len <= n) {
            ret = add_mod(ret, cal(nodes[r].s, nodes[r].e, idx), mod);
        }
        idx += len;

        if (nodes[r].child[1]) ret = add_mod(ret, solve_impl(nodes[r].child[1]), mod);
        return ret;
    }

    i64 solve() {
        idx = -1;
        return solve_impl(nodes[0].child[0]);
    }

    void split_node(int id, i64 offset) {
        if (offset <= 1) return;
        splay(id, 0);
        int dir = nodes[id].s <= nodes[id].e ? 1 : -1;
        i64 olde = nodes[id].e;
        nodes[id].e = nodes[id].s + dir * (offset - 2);

        int z = nodes[id].child[1];
        push_down(z);
        while (nodes[z].child[0]) {
            z = nodes[z].child[0];
            push_down(z);
        }

        int new_id = top++;
        nodes[new_id].parent = z;
        nodes[new_id].child[0] = nodes[new_id].child[1] = 0;
        nodes[new_id].s = nodes[id].e + dir;
        nodes[new_id].e = olde;
        nodes[new_id].cnt = seg_len(nodes[new_id]);
        nodes[new_id].rev = false;
        nodes[z].child[0] = new_id;

        while (z != 0) {
            update(z);
            z = nodes[z].parent;
        }
    }

    void reverse_range(i64 s, i64 e) {
        if (s > e) std::swap(s, e);

        for (i64 target : {s + 2, e + 3}) {
            i64 x = target;
            int id = find(x);
            if (x != 1) {
                split_node(id, x);
            }
        }

        i64 x = s + 1;
        const int id1 = find(x);
        splay(id1, 0);

        i64 y = e + 3;
        const int id2 = find(y);
        splay(id2, id1);

        const int id3 = nodes[id2].child[0];
        nodes[id3].rev = !nodes[id3].rev;
        std::swap(nodes[id3].s, nodes[id3].e);
    }
};

u128 brute_force(i64 n, i64 k) {
    std::vector<i64> arr(static_cast<std::size_t>(n));
    for (i64 i = 0; i < n; ++i) arr[static_cast<std::size_t>(i)] = i;
    i64 s = 1 % n;
    i64 t = 1 % n;
    for (i64 step = 0; step < k; ++step) {
        i64 l = std::min(s, t);
        i64 r = std::max(s, t);
        std::reverse(arr.begin() + l, arr.begin() + r + 1);
        i64 ns = add_mod_n(s, t, n);
        i64 nt = add_mod_n(t, ns, n);
        s = ns;
        t = nt;
    }
    u128 sum = 0;
    for (i64 i = 0; i < n; ++i) sum += static_cast<u128>(i) * static_cast<u128>(arr[static_cast<std::size_t>(i)]);
    return sum;
}

i64 compute_fast(i64 n, i64 k, i64 mod) {
    SplayTree st(n, mod);
    i64 s = 1 % n;
    i64 t = 1 % n;
    for (i64 step = 0; step < k; ++step) {
        if (s != t) {
            st.reverse_range(s, t);
        }
        const i64 next_s = add_mod_n(s, t, n);
        const i64 next_t = add_mod_n(t, next_s, n);
        s = next_s;
        t = next_t;
    }
    return st.solve();
}

void run_checkpoints(i64 mod) {
    const u128 v1 = brute_force(5, 4);
    const u128 v2 = brute_force(100, 100);
    const u128 v3 = brute_force(10000, 10000);

    if (v1 != static_cast<u128>(27)) throw std::runtime_error("checkpoint R(5,4)");
    if (v2 != static_cast<u128>(246597)) throw std::runtime_error("checkpoint R(100,100)");
    if (v3 != static_cast<u128>(249275481640LL)) throw std::runtime_error("checkpoint R(10000,10000)");

    if (compute_fast(5, 4, mod) != static_cast<i64>(v1 % mod)) {
        throw std::runtime_error("fast mismatch R(5,4)");
    }
    if (compute_fast(100, 100, mod) != static_cast<i64>(v2 % mod)) {
        throw std::runtime_error("fast mismatch R(100,100)");
    }
    if (compute_fast(10000, 10000, mod) != static_cast<i64>(v3 % mod)) {
        throw std::runtime_error("fast mismatch R(10000,10000)");
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) return 1;
        if (options.run_checkpoints) {
            run_checkpoints(options.mod);
        }

        std::cout << compute_fast(options.n, options.k, options.mod) << '\n';
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
