#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct ShapeInfo {
    int left = -1;
    int right = -1;
    int leaves = 1;
};

struct HeapItem {
    u64 prod;
    int i;
    int j;
};

struct HeapCmp {
    bool operator()(const HeapItem& a, const HeapItem& b) const {
        if (a.prod != b.prod) return a.prod > b.prod;
        if (a.i != b.i) return a.i > b.i;
        return a.j > b.j;
    }
};

struct Sequence {
    bool needed = false;
    bool leaf = false;
    bool started = false;
    bool exhausted = false;

    int left = -1;
    int right = -1;
    int next_i = 0;

    u64 next_prime = 2;

    std::vector<u64> values;

    std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
    std::unordered_set<u64> in_heap;
};

static std::mt19937_64 rng(0);

static std::vector<ShapeInfo> shapes(1);  // id 0 = LEAF
static std::map<std::pair<int, int>, int> shape_id_by_children;
static std::unordered_map<u64, int> shape_cache;
static std::unordered_map<u64, u64> best_divisor_cache;

static std::unordered_map<u64, bool> prime_cache;
static std::unordered_map<u64, std::map<u64, int>> factor_cache;

static u64 MAXVAL = 0;

static std::vector<Sequence> seqs;

static u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) r = mul_mod(r, a, mod);
        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool is_prime64(u64 n) {
    auto it = prime_cache.find(n);
    if (it != prime_cache.end()) return it->second;

    if (n < 2) return prime_cache[n] = false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return prime_cache[n] = (n == p);
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return prime_cache[n] = false;
    }

    return prime_cache[n] = true;
}

static u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    if (n % 3ULL == 0) return 3;

    std::uniform_int_distribution<u64> dist(2, n - 2);

    while (true) {
        u64 c = dist(rng);
        u64 x = dist(rng);
        u64 y = x;
        u64 d = 1;

        while (d == 1) {
            x = (mul_mod(x, x, n) + c) % n;
            y = (mul_mod(y, y, n) + c) % n;
            y = (mul_mod(y, y, n) + c) % n;
            u64 diff = (x > y ? x - y : y - x);
            d = std::gcd(diff, n);
        }

        if (d != n) return d;
    }
}

static void factor_rec(u64 n, std::map<u64, int>& fac) {
    if (n == 1) return;
    if (is_prime64(n)) {
        fac[n]++;
        return;
    }
    u64 d = pollard_rho(n);
    factor_rec(d, fac);
    factor_rec(n / d, fac);
}

static std::map<u64, int> factorize(u64 n) {
    auto it = factor_cache.find(n);
    if (it != factor_cache.end()) return it->second;
    std::map<u64, int> fac;
    factor_rec(n, fac);
    factor_cache[n] = fac;
    return fac;
}

static void gen_divisors_rec(const std::vector<std::pair<u64, int>>& f, int idx, u64 cur, std::vector<u64>& out) {
    if (idx == static_cast<int>(f.size())) {
        out.push_back(cur);
        return;
    }
    auto [p, e] = f[idx];
    u64 v = 1;
    for (int i = 0; i <= e; ++i) {
        gen_divisors_rec(f, idx + 1, cur * v, out);
        v *= p;
    }
}

static u64 best_divisor_le_sqrt(u64 n) {
    auto it = best_divisor_cache.find(n);
    if (it != best_divisor_cache.end()) return it->second;

    auto fac_map = factorize(n);
    std::vector<std::pair<u64, int>> f(fac_map.begin(), fac_map.end());
    std::vector<u64> divs;
    divs.reserve(1024);
    gen_divisors_rec(f, 0, 1, divs);

    u64 root = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((root + 1) * (root + 1) <= n) ++root;
    while (root * root > n) --root;

    u64 best = 1;
    for (u64 d : divs) {
        if (d <= root && d > best) best = d;
    }

    best_divisor_cache[n] = best;
    return best;
}

static int make_shape(int l, int r) {
    auto key = std::make_pair(l, r);
    auto it = shape_id_by_children.find(key);
    if (it != shape_id_by_children.end()) return it->second;

    ShapeInfo s;
    s.left = l;
    s.right = r;
    s.leaves = shapes[l].leaves + shapes[r].leaves;
    int id = static_cast<int>(shapes.size());
    shapes.push_back(s);
    shape_id_by_children[key] = id;
    return id;
}

static int shape_of(u64 n) {
    auto it = shape_cache.find(n);
    if (it != shape_cache.end()) return it->second;
    if (is_prime64(n)) {
        shape_cache[n] = 0;
        return 0;
    }

    u64 d = best_divisor_le_sqrt(n);
    u64 a = d;
    u64 b = n / d;
    if (a > b) std::swap(a, b);

    int l = shape_of(a);
    int r = shape_of(b);
    int id = make_shape(l, r);
    shape_cache[n] = id;
    return id;
}

static std::optional<u64> next_value(int sh);

static std::optional<u64> get_value_maybe(int sh, int idx) {
    auto& seq = seqs[sh];
    while (static_cast<int>(seq.values.size()) <= idx) {
        auto nv = next_value(sh);
        if (!nv.has_value()) return std::nullopt;
        seq.values.push_back(*nv);
    }
    return seq.values[idx];
}

static int ensure_right_ge(int right_sh, u64 x) {
    auto& rseq = seqs[right_sh];
    if (rseq.exhausted && (rseq.values.empty() || rseq.values.back() < x)) return -1;

    while (rseq.values.empty() || rseq.values.back() < x) {
        auto nv = next_value(right_sh);
        if (!nv.has_value()) {
            rseq.exhausted = true;
            break;
        }
        rseq.values.push_back(*nv);
    }

    if (rseq.values.empty() || rseq.values.back() < x) return -1;
    auto it = std::lower_bound(rseq.values.begin(), rseq.values.end(), x);
    return static_cast<int>(it - rseq.values.begin());
}

static void push_pair(int sh, int i, int j) {
    auto& seq = seqs[sh];
    u64 key = (static_cast<u64>(i) << 32) | static_cast<u64>(j);
    if (seq.in_heap.count(key)) return;

    auto x = get_value_maybe(seq.left, i);
    auto y = get_value_maybe(seq.right, j);
    if (!x.has_value() || !y.has_value()) return;
    if (*x > *y) return;

    u128 prod = static_cast<u128>(*x) * static_cast<u128>(*y);
    if (prod > MAXVAL) return;

    seq.in_heap.insert(key);
    seq.heap.push({static_cast<u64>(prod), i, j});
}

static void start_node(int sh) {
    auto& seq = seqs[sh];
    if (seq.started) return;
    seq.started = true;

    auto x0 = get_value_maybe(seq.left, 0);
    if (!x0.has_value()) {
        seq.exhausted = true;
        return;
    }
    int j0 = ensure_right_ge(seq.right, *x0);
    if (j0 != -1) {
        push_pair(sh, 0, j0);
    }
    seq.next_i = 1;
}

static void maybe_add_more_i(int sh, u64 current_min) {
    auto& seq = seqs[sh];
    while (true) {
        auto x = get_value_maybe(seq.left, seq.next_i);
        if (!x.has_value()) return;

        if (!seq.heap.empty() && static_cast<u128>(*x) * static_cast<u128>(*x) > current_min) {
            return;
        }

        int j0 = ensure_right_ge(seq.right, *x);
        if (j0 != -1) {
            push_pair(sh, seq.next_i, j0);
        }
        ++seq.next_i;
    }
}

static std::optional<u64> next_candidate_product(int sh) {
    auto& seq = seqs[sh];
    start_node(sh);
    if (seq.exhausted) return std::nullopt;

    while (true) {
        if (seq.heap.empty()) {
            auto x = get_value_maybe(seq.left, seq.next_i);
            if (!x.has_value()) {
                seq.exhausted = true;
                return std::nullopt;
            }
            int j0 = ensure_right_ge(seq.right, *x);
            if (j0 == -1) {
                seq.exhausted = true;
                return std::nullopt;
            }
            push_pair(sh, seq.next_i, j0);
            ++seq.next_i;
            continue;
        }

        HeapItem cur = seq.heap.top();
        seq.heap.pop();
        u64 key = (static_cast<u64>(cur.i) << 32) | static_cast<u64>(cur.j);
        seq.in_heap.erase(key);

        maybe_add_more_i(sh, cur.prod);
        push_pair(sh, cur.i, cur.j + 1);

        return cur.prod;
    }
}

static std::optional<u64> next_value(int sh) {
    auto& seq = seqs[sh];
    if (seq.exhausted) return std::nullopt;

    if (seq.leaf) {
        if (seq.next_prime == 2) {
            seq.next_prime = 3;
            return 2;
        }
        for (u64 p = seq.next_prime; p <= MAXVAL; p += 2) {
            if (is_prime64(p)) {
                seq.next_prime = p + 2;
                return p;
            }
        }
        seq.exhausted = true;
        return std::nullopt;
    }

    while (true) {
        auto cand = next_candidate_product(sh);
        if (!cand.has_value()) {
            seq.exhausted = true;
            return std::nullopt;
        }
        if (!seq.values.empty() && *cand == seq.values.back()) {
            continue;
        }
        if (shape_of(*cand) == sh) {
            return cand;
        }
    }
}

static u64 get_value(int sh, int idx) {
    auto v = get_value_maybe(sh, idx);
    if (!v.has_value()) {
        throw std::runtime_error("sequence exhausted");
    }
    return *v;
}

static u64 double_factorial(int n) {
    u64 r = 1;
    for (int k = n; k >= 2; k -= 2) {
        r *= static_cast<u64>(k);
    }
    return r;
}

int main() {
    const int max_n = 31;
    MAXVAL = double_factorial(max_n);

    std::vector<int> targets;
    targets.reserve(max_n - 1);
    for (int n = 2; n <= max_n; ++n) {
        targets.push_back(shape_of(double_factorial(n)));
    }

    std::vector<bool> needed(shapes.size(), false);
    std::function<void(int)> collect = [&](int sh) {
        if (needed[sh]) return;
        needed[sh] = true;
        if (sh == 0) return;
        collect(shapes[sh].left);
        collect(shapes[sh].right);
    };
    for (int sh : targets) collect(sh);

    seqs.resize(shapes.size());
    for (int sh = 0; sh < static_cast<int>(shapes.size()); ++sh) {
        if (!needed[sh]) continue;
        seqs[sh].needed = true;
        seqs[sh].leaf = (sh == 0);
        if (sh != 0) {
            seqs[sh].left = shapes[sh].left;
            seqs[sh].right = shapes[sh].right;
        }
    }

    std::vector<u64> M(max_n + 1, 0);
    for (int n = 2; n <= max_n; ++n) {
        M[n] = get_value(targets[n - 2], 0);
    }

    assert(M[9] == 72ULL);

    u128 sum = 0;
    for (int n = 2; n <= max_n; ++n) {
        sum += M[n];
    }

    std::cout << static_cast<u64>(sum) << '\n';
    return 0;
}
