#include <cassert>
#include <cstdint>
#include <future>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

namespace {
constexpr uint64_t MOD = 1000000000ULL;

const vector<vector<string>> kBaseFactors = {
    {},
    {"0", "1"},
    {"00", "01", "10", "11"},
    {"001", "010", "011", "100", "101", "110"},
};

thread_local unordered_map<long long, unsigned long long> p_cache;
thread_local unordered_map<long long, unsigned long long> s_cache;

// For n>=4: p(2n)=p(n)+p(n+1), p(2n+1)=2*p(n+1).
// For n>=4: S(2n)=3*S(n)+S(n+1)-8, S(2n+1)=4*S(n)+3*p(n+1)-8.
unsigned long long p_count(long long n) {
    if (n <= 0) return 0;
    if (n == 1) return 2;
    if (n == 2) return 4;
    if (n == 3) return 6;
    auto it = p_cache.find(n);
    if (it != p_cache.end()) return it->second;
    unsigned long long res = 0;
    if ((n & 1LL) == 0) {
        long long m = n / 2;
        res = p_count(m) + p_count(m + 1);
    } else {
        long long m = n / 2;
        res = 2ULL * p_count(m + 1);
    }
    p_cache[n] = res;
    return res;
}

unsigned long long s_count(long long n) {
    if (n <= 0) return 0;
    if (n == 1) return 2;
    if (n == 2) return 6;
    if (n == 3) return 12;
    auto it = s_cache.find(n);
    if (it != s_cache.end()) return it->second;
    unsigned long long res = 0;
    if ((n & 1LL) == 0) {
        long long m = n / 2;
        __int128 tmp = 3; 
        tmp = tmp * s_count(m) + s_count(m + 1) - 8;
        res = static_cast<unsigned long long>(tmp);
    } else {
        long long m = n / 2;
        __int128 tmp = 4;
        tmp = tmp * s_count(m) + 3 * static_cast<__int128>(p_count(m + 1)) - 8;
        res = static_cast<unsigned long long>(tmp);
    }
    s_cache[n] = res;
    return res;
}

unsigned long long count_upto(long long len) {
    if (len <= 0) return 0;
    return s_count(len) / 2ULL;
}

long long find_length_for_index(unsigned long long n) {
    long long lo = 1;
    long long hi = 1;
    while (count_upto(hi) < n) {
        hi *= 2;
    }
    while (lo < hi) {
        long long mid = lo + (hi - lo) / 2;
        if (count_upto(mid) >= n) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return lo;
}

struct WordBuilder {
    struct Node {
        enum Type { EMPTY, LEAF, CONCAT, MU } type = EMPTY;
        int left = -1;
        int right = -1;
        int child = -1;
        string bits;
        long long len = 0;
        int first = -1;
        int last = -1;
    };

    uint64_t mod = MOD;
    int max_i = 64;
    vector<uint64_t> base;
    vector<unordered_map<long long, pair<uint64_t, uint64_t>>> pow_sum_cache;

    vector<Node> nodes;
    vector<vector<uint64_t>> eval_cache;

    unordered_map<uint64_t, int> prefix_cache;
    unordered_map<uint64_t, int> suffix_cache;
    unordered_map<uint64_t, int> middle_cache;

    int empty_id = 0;
    int leaf0_id = 0;
    int leaf1_id = 0;

    explicit WordBuilder(uint64_t mod_, int max_i_ = 64) : mod(mod_), max_i(max_i_) {
        base.assign(max_i + 2, 0);
        base[0] = 2 % mod;
        for (int i = 1; i < static_cast<int>(base.size()); ++i) {
            base[i] = mul_mod(base[i - 1], base[i - 1]);
        }
        pow_sum_cache.resize(base.size());

        nodes.reserve(512);
        eval_cache.reserve(512);
        Node empty;
        empty.type = Node::EMPTY;
        empty.len = 0;
        nodes.push_back(empty);
        eval_cache.emplace_back();
        empty_id = 0;

        Node leaf0;
        leaf0.type = Node::LEAF;
        leaf0.bits = "0";
        leaf0.len = 1;
        leaf0.first = 0;
        leaf0.last = 0;
        nodes.push_back(leaf0);
        eval_cache.emplace_back();
        leaf0_id = static_cast<int>(nodes.size() - 1);

        Node leaf1;
        leaf1.type = Node::LEAF;
        leaf1.bits = "1";
        leaf1.len = 1;
        leaf1.first = 1;
        leaf1.last = 1;
        nodes.push_back(leaf1);
        eval_cache.emplace_back();
        leaf1_id = static_cast<int>(nodes.size() - 1);
    }

    uint64_t mul_mod(uint64_t a, uint64_t b) const {
        return static_cast<uint64_t>((__int128)a * b % mod);
    }

    int make_leaf(const string& bits) {
        if (bits.empty()) return empty_id;
        if (bits.size() == 1) return bits[0] == '0' ? leaf0_id : leaf1_id;
        Node node;
        node.type = Node::LEAF;
        node.bits = bits;
        node.len = static_cast<long long>(bits.size());
        node.first = bits.front() - '0';
        node.last = bits.back() - '0';
        nodes.push_back(node);
        eval_cache.emplace_back();
        return static_cast<int>(nodes.size() - 1);
    }

    int make_bit(int b) { return b == 0 ? leaf0_id : leaf1_id; }

    int make_concat(int left, int right) {
        if (nodes[left].len == 0) return right;
        if (nodes[right].len == 0) return left;
        Node node;
        node.type = Node::CONCAT;
        node.left = left;
        node.right = right;
        node.len = nodes[left].len + nodes[right].len;
        node.first = nodes[left].len ? nodes[left].first : nodes[right].first;
        node.last = nodes[right].len ? nodes[right].last : nodes[left].last;
        nodes.push_back(node);
        eval_cache.emplace_back();
        return static_cast<int>(nodes.size() - 1);
    }

    int make_mu(int child) {
        if (nodes[child].len == 0) return empty_id;
        Node node;
        node.type = Node::MU;
        node.child = child;
        node.len = nodes[child].len * 2;
        node.first = nodes[child].first;
        node.last = 1 - nodes[child].last;
        nodes.push_back(node);
        eval_cache.emplace_back();
        return static_cast<int>(nodes.size() - 1);
    }

    int prefix(int node_id) {
        const Node& node = nodes[node_id];
        if (node.len <= 1) return empty_id;
        uint64_t key = (static_cast<uint64_t>(node_id) << 2) | 0ULL;
        auto it = prefix_cache.find(key);
        if (it != prefix_cache.end()) return it->second;
        int res = empty_id;
        if (node.type == Node::LEAF) {
            res = make_leaf(node.bits.substr(0, static_cast<size_t>(node.len - 1)));
        } else if (node.type == Node::CONCAT) {
            const Node& right = nodes[node.right];
            if (right.len > 1) {
                res = make_concat(node.left, prefix(node.right));
            } else if (right.len == 1) {
                res = node.left;
            } else {
                res = prefix(node.left);
            }
        } else if (node.type == Node::MU) {
            int p = prefix(node.child);
            int mu_p = make_mu(p);
            int bit = make_bit(nodes[node.child].last);
            res = make_concat(mu_p, bit);
        }
        prefix_cache[key] = res;
        return res;
    }

    int suffix(int node_id) {
        const Node& node = nodes[node_id];
        if (node.len <= 1) return empty_id;
        uint64_t key = (static_cast<uint64_t>(node_id) << 2) | 1ULL;
        auto it = suffix_cache.find(key);
        if (it != suffix_cache.end()) return it->second;
        int res = empty_id;
        if (node.type == Node::LEAF) {
            res = make_leaf(node.bits.substr(1));
        } else if (node.type == Node::CONCAT) {
            const Node& left = nodes[node.left];
            if (left.len > 1) {
                res = make_concat(suffix(node.left), node.right);
            } else if (left.len == 1) {
                res = node.right;
            } else {
                res = suffix(node.right);
            }
        } else if (node.type == Node::MU) {
            int s = suffix(node.child);
            int mu_s = make_mu(s);
            int bit = make_bit(1 - nodes[node.child].first);
            res = make_concat(bit, mu_s);
        }
        suffix_cache[key] = res;
        return res;
    }

    int middle(int node_id) {
        const Node& node = nodes[node_id];
        if (node.len <= 2) return empty_id;
        uint64_t key = (static_cast<uint64_t>(node_id) << 2) | 2ULL;
        auto it = middle_cache.find(key);
        if (it != middle_cache.end()) return it->second;
        int res = empty_id;
        if (node.type == Node::LEAF) {
            res = make_leaf(node.bits.substr(1, static_cast<size_t>(node.len - 2)));
        } else if (node.type == Node::CONCAT) {
            int tmp = suffix(node_id);
            res = prefix(tmp);
        } else if (node.type == Node::MU) {
            int m = middle(node.child);
            int mu_m = make_mu(m);
            int bit_left = make_bit(1 - nodes[node.child].first);
            int bit_right = make_bit(nodes[node.child].last);
            res = make_concat(bit_left, make_concat(mu_m, bit_right));
        }
        middle_cache[key] = res;
        return res;
    }

    int odd_even_map(int ancestor) {
        int mid = middle(ancestor);
        int mu_mid = make_mu(mid);
        int bit_left = make_bit(1 - nodes[ancestor].first);
        int bit_right = make_bit(nodes[ancestor].last);
        return make_concat(bit_left, make_concat(mu_mid, bit_right));
    }

    int odd_odd_map(int ancestor) {
        int suf = suffix(ancestor);
        int mu_suf = make_mu(suf);
        int bit_left = make_bit(1 - nodes[ancestor].first);
        return make_concat(bit_left, mu_suf);
    }

    int even_odd_map(int ancestor) {
        int pref = prefix(ancestor);
        int mu_pref = make_mu(pref);
        int bit_right = make_bit(nodes[ancestor].last);
        return make_concat(mu_pref, bit_right);
    }

    pair<uint64_t, uint64_t> pow_sum(int i, long long len) {
        auto& cache = pow_sum_cache[i];
        auto it = cache.find(len);
        if (it != cache.end()) return it->second;
        pair<uint64_t, uint64_t> res;
        if (len == 0) {
            res = {1 % mod, 0};
        } else if (len == 1) {
            res = {base[i] % mod, 1 % mod};
        } else if ((len & 1LL) == 0) {
            auto half = pow_sum(i, len / 2);
            uint64_t p = half.first;
            uint64_t s = half.second;
            uint64_t p2 = mul_mod(p, p);
            uint64_t s2 = mul_mod(s, (p + 1) % mod);
            res = {p2, s2};
        } else {
            auto prev = pow_sum(i, len - 1);
            uint64_t p = prev.first;
            uint64_t s = prev.second;
            uint64_t p2 = mul_mod(p, base[i]);
            uint64_t s2 = (s + p) % mod;
            res = {p2, s2};
        }
        cache[len] = res;
        return res;
    }

    uint64_t pow_base(int i, long long len) {
        return pow_sum(i, len).first;
    }

    uint64_t sum_geom(int i, long long len) {
        return pow_sum(i, len).second;
    }

    const vector<uint64_t>& eval(int node_id) {
        auto& cached = eval_cache[node_id];
        if (!cached.empty()) return cached;
        cached.assign(max_i + 2, 0);
        const Node& node = nodes[node_id];
        if (node.len == 0) {
            return cached;
        }
        if (node.type == Node::LEAF) {
            for (int i = 0; i < static_cast<int>(cached.size()); ++i) {
                uint64_t v = 0;
                for (char c : node.bits) {
                    v = (mul_mod(v, base[i]) + static_cast<uint64_t>(c - '0')) % mod;
                }
                cached[i] = v;
            }
        } else if (node.type == Node::CONCAT) {
            const auto& left = eval(node.left);
            const auto& right = eval(node.right);
            long long len_r = nodes[node.right].len;
            for (int i = 0; i < static_cast<int>(cached.size()); ++i) {
                uint64_t term = mul_mod(left[i], pow_base(i, len_r));
                cached[i] = (term + right[i]) % mod;
            }
        } else if (node.type == Node::MU) {
            const auto& child = eval(node.child);
            long long len_c = nodes[node.child].len;
            for (int i = 0; i + 1 < static_cast<int>(cached.size()); ++i) {
                uint64_t sum = sum_geom(i + 1, len_c);
                uint64_t factor = (base[i] + mod - 1) % mod;
                uint64_t term = mul_mod(factor, child[i + 1]);
                cached[i] = (sum + term) % mod;
            }
        }
        return cached;
    }

    string materialize(int node_id, size_t limit) {
        const Node& node = nodes[node_id];
        if (node.len > static_cast<long long>(limit)) return string();
        if (node.type == Node::EMPTY) return string();
        if (node.type == Node::LEAF) return node.bits;
        if (node.type == Node::CONCAT) {
            string left = materialize(node.left, limit);
            if (left.empty() && nodes[node.left].len > 0) return string();
            string right = materialize(node.right, limit);
            if (right.empty() && nodes[node.right].len > 0) return string();
            return left + right;
        }
        string child = materialize(node.child, limit);
        if (child.empty() && nodes[node.child].len > 0) return string();
        string out;
        out.reserve(child.size() * 2);
        for (char c : child) {
            if (c == '0') out += "01";
            else out += "10";
        }
        return out;
    }
};

int select_node(long long len, unsigned long long k, WordBuilder& builder) {
    if (len <= 3) {
        return builder.make_leaf(kBaseFactors[static_cast<size_t>(len)][static_cast<size_t>(k - 1)]);
    }
    // Lexicographic split: O1 (odd-aligned from leading 1 ancestors), E (even-aligned), O0 (odd-aligned from leading 0 ancestors).
    if ((len & 1LL) == 0) {
        long long n = len / 2;
        unsigned long long size_O1 = p_count(n + 1) / 2ULL;
        unsigned long long size_E = p_count(n);
        if (k <= size_O1) {
            int ancestor = select_node(n + 1, p_count(n + 1) / 2ULL + k, builder);
            return builder.odd_even_map(ancestor);
        }
        if (k <= size_O1 + size_E) {
            int ancestor = select_node(n, k - size_O1, builder);
            return builder.make_mu(ancestor);
        }
        int ancestor = select_node(n + 1, k - size_O1 - size_E, builder);
        return builder.odd_even_map(ancestor);
    }
    long long n = len / 2;
    unsigned long long size_O1 = p_count(n + 1) / 2ULL;
    unsigned long long size_E = p_count(n + 1);
    if (k <= size_O1) {
        int ancestor = select_node(n + 1, p_count(n + 1) / 2ULL + k, builder);
        return builder.odd_odd_map(ancestor);
    }
    if (k <= size_O1 + size_E) {
        int ancestor = select_node(n + 1, k - size_O1, builder);
        return builder.even_odd_map(ancestor);
    }
    int ancestor = select_node(n + 1, k - size_O1 - size_E, builder);
    return builder.odd_odd_map(ancestor);
}

uint64_t compute_A_mod(unsigned long long n) {
    if (n == 0) return 0;
    long long len = find_length_for_index(n);
    unsigned long long prev = count_upto(len - 1);
    unsigned long long k = n - prev;
    unsigned long long offset = p_count(len) / 2ULL + k;
    WordBuilder builder(MOD);
    int root = select_node(len, offset, builder);
    const auto& vals = builder.eval(root);
    return vals[0] % MOD;
}

void run_validation() {
    assert(p_count(1) == 2);
    assert(p_count(2) == 4);
    assert(p_count(3) == 6);
    assert(p_count(4) == 10);
    assert(s_count(3) == 12);
    assert(s_count(4) == 22);

    {
        unsigned long long n = 100;
        long long len = find_length_for_index(n);
        unsigned long long prev = count_upto(len - 1);
        unsigned long long k = n - prev;
        unsigned long long offset = p_count(len) / 2ULL + k;
        WordBuilder builder(MOD);
        int root = select_node(len, offset, builder);
        string bits = builder.materialize(root, 128);
        assert(!bits.empty());
        unsigned long long value = 0;
        for (char c : bits) value = (value << 1) + (c - '0');
        assert(value == 3251ULL);
    }
    {
        unsigned long long n = 1000;
        long long len = find_length_for_index(n);
        unsigned long long prev = count_upto(len - 1);
        unsigned long long k = n - prev;
        unsigned long long offset = p_count(len) / 2ULL + k;
        WordBuilder builder(MOD);
        int root = select_node(len, offset, builder);
        string bits = builder.materialize(root, 256);
        assert(!bits.empty());
        unsigned long long value = 0;
        for (char c : bits) value = (value << 1) + (c - '0');
        assert(value == 80852364498ULL);
    }
}

} // namespace

int main() {
    run_validation();

    vector<unsigned long long> targets;
    unsigned long long v = 1;
    for (int k = 1; k <= 18; ++k) {
        v *= 10ULL;
        targets.push_back(v);
    }

    vector<future<uint64_t>> futures;
    futures.reserve(targets.size());
    for (unsigned long long n : targets) {
        futures.emplace_back(async(launch::async, [n]() { return compute_A_mod(n); }));
    }

    uint64_t sum_mod = 0;
    for (auto& fut : futures) {
        sum_mod = (sum_mod + fut.get()) % MOD;
    }

    cout << setw(9) << setfill('0') << sum_mod << "\n";
    return 0;
}
