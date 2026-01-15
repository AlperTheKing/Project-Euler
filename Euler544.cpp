#include <array>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <vector>

using namespace std;

constexpr int64_t MOD = 1000000007LL;
constexpr int MAX_DEG = 90;
constexpr int MAX_FRONTIER = 10;

using Poly = array<int64_t, MAX_DEG + 1>;

struct State {
    uint8_t size = 0;
    uint8_t num_labels = 0;
    array<uint8_t, MAX_FRONTIER + 1> labels{};
};

static uint64_t encode_state(const State& s) {
    uint64_t key = s.size;
    for (uint8_t i = 0; i < s.size; ++i) {
        key = (key << 4) | s.labels[i];
    }
    return key;
}

static State canonicalize(const State& s) {
    State out = s;
    array<uint8_t, 16> remap;
    remap.fill(0xFF);
    uint8_t next = 0;
    for (uint8_t i = 0; i < out.size; ++i) {
        uint8_t lbl = out.labels[i];
        if (remap[lbl] == 0xFF) remap[lbl] = next++;
        out.labels[i] = remap[lbl];
    }
    out.num_labels = next;
    return out;
}

static inline void poly_add(Poly& dest, const Poly& src) {
    for (int i = 0; i <= MAX_DEG; ++i) {
        int64_t v = dest[i] + src[i];
        if (v >= MOD) v -= MOD;
        dest[i] = v;
    }
}

static inline void poly_sub(Poly& dest, const Poly& src) {
    for (int i = 0; i <= MAX_DEG; ++i) {
        int64_t v = dest[i] - src[i];
        if (v < 0) v += MOD;
        dest[i] = v;
    }
}

static inline void poly_mul_q(Poly& p) {
    for (int i = MAX_DEG - 1; i >= 0; --i) {
        p[i + 1] = p[i];
    }
    p[0] = 0;
}

struct StateMap {
    vector<State> states;
    vector<Poly> polys;
    unordered_map<uint64_t, size_t> index;

    void reserve(size_t n) {
        states.reserve(n);
        polys.reserve(n);
        index.reserve(n * 2);
    }

    void add_state(const State& s, const Poly& p, bool negate) {
        uint64_t key = encode_state(s);
        auto it = index.find(key);
        if (it == index.end()) {
            states.push_back(s);
            Poly poly = p;
            if (negate) {
                for (int i = 0; i <= MAX_DEG; ++i) {
                    if (poly[i] != 0) poly[i] = MOD - poly[i];
                }
            }
            polys.push_back(poly);
            index.emplace(key, states.size() - 1);
        } else {
            Poly& dest = polys[it->second];
            if (negate) {
                poly_sub(dest, p);
            } else {
                poly_add(dest, p);
            }
        }
    }

    void merge_from(const StateMap& other) {
        for (size_t i = 0; i < other.states.size(); ++i) {
            add_state(other.states[i], other.polys[i], false);
        }
    }
};

constexpr size_t kParallelThreshold = 2000;

template <typename Func>
static StateMap transform_parallel(const StateMap& cur, int threads, size_t reserve_hint, Func func) {
    if (threads <= 1 || cur.states.size() < kParallelThreshold) {
        StateMap out;
        out.reserve(reserve_hint);
        for (size_t i = 0; i < cur.states.size(); ++i) {
            func(cur.states[i], cur.polys[i], out);
        }
        return out;
    }

    int t = min<int>(threads, static_cast<int>(cur.states.size()));
    vector<StateMap> locals(static_cast<size_t>(t));
    vector<thread> workers;
    workers.reserve(static_cast<size_t>(t));

    size_t total = cur.states.size();
    size_t chunk = (total + static_cast<size_t>(t) - 1) / static_cast<size_t>(t);

    for (int ti = 0; ti < t; ++ti) {
        size_t start = static_cast<size_t>(ti) * chunk;
        size_t end = min(total, start + chunk);
        if (start >= end) continue;
        locals[static_cast<size_t>(ti)].reserve((end - start) * 2);
        workers.emplace_back([&, start, end, ti]() {
            StateMap& out = locals[static_cast<size_t>(ti)];
            for (size_t i = start; i < end; ++i) {
                func(cur.states[i], cur.polys[i], out);
            }
        });
    }

    for (auto& w : workers) w.join();

    StateMap merged;
    merged.reserve(reserve_hint);
    for (const auto& local : locals) {
        merged.merge_from(local);
    }
    return merged;
}

static StateMap add_vertex_end(const StateMap& cur) {
    StateMap out;
    out.reserve(cur.states.size());
    for (size_t i = 0; i < cur.states.size(); ++i) {
        State ns = cur.states[i];
        ns.labels[ns.size] = ns.num_labels;
        ns.size++;
        ns.num_labels++;
        out.add_state(ns, cur.polys[i], false);
    }
    return out;
}

static StateMap add_edge(const StateMap& cur, int pos_a, int pos_b, int threads) {
    size_t reserve_hint = cur.states.size() * 2;
    auto func = [pos_a, pos_b](const State& st, const Poly& poly, StateMap& out) {
        out.add_state(st, poly, false);
        State merged = st;
        uint8_t la = merged.labels[static_cast<size_t>(pos_a)];
        uint8_t lb = merged.labels[static_cast<size_t>(pos_b)];
        if (la != lb) {
            for (uint8_t i = 0; i < merged.size; ++i) {
                if (merged.labels[i] == lb) merged.labels[i] = la;
            }
            merged = canonicalize(merged);
        }
        out.add_state(merged, poly, true);
    };
    return transform_parallel(cur, threads, reserve_hint, func);
}

static StateMap remove_pos_swap(const StateMap& cur, int pos, int threads) {
    size_t reserve_hint = cur.states.size();
    auto func = [pos](const State& st, const Poly& poly, StateMap& out) {
        uint8_t label = st.labels[static_cast<size_t>(pos)];
        bool unique = true;
        for (uint8_t i = 0; i < st.size; ++i) {
            if (i != pos && st.labels[i] == label) {
                unique = false;
                break;
            }
        }
        State ns = st;
        uint8_t last = static_cast<uint8_t>(ns.size - 1);
        ns.labels[static_cast<size_t>(pos)] = ns.labels[last];
        ns.size--;
        ns = canonicalize(ns);
        Poly p = poly;
        if (unique) poly_mul_q(p);
        out.add_state(ns, p, false);
    };
    return transform_parallel(cur, threads, reserve_hint, func);
}

static StateMap remove_last(const StateMap& cur, int threads) {
    size_t reserve_hint = cur.states.size();
    auto func = [](const State& st, const Poly& poly, StateMap& out) {
        uint8_t pos = static_cast<uint8_t>(st.size - 1);
        uint8_t label = st.labels[pos];
        bool unique = true;
        for (uint8_t i = 0; i + 1 < st.size; ++i) {
            if (st.labels[i] == label) {
                unique = false;
                break;
            }
        }
        State ns = st;
        ns.size--;
        ns = canonicalize(ns);
        Poly p = poly;
        if (unique) poly_mul_q(p);
        out.add_state(ns, p, false);
    };
    return transform_parallel(cur, threads, reserve_hint, func);
}

static Poly chromatic_poly_grid(int rows, int cols, int threads) {
    if (rows > cols) swap(rows, cols);
    StateMap cur;
    cur.reserve(1);

    State empty;
    Poly base{};
    base[0] = 1;
    cur.add_state(empty, base, false);

    int frontier = 0;
    for (int c = 0; c < cols; ++c) {
        for (int r = 0; r < rows; ++r) {
            cur = add_vertex_end(cur);
            frontier += 1;
            int new_pos = frontier - 1;
            if (c > 0) {
                cur = add_edge(cur, r, new_pos, threads);
            }
            if (r > 0) {
                cur = add_edge(cur, r - 1, new_pos, threads);
            }
            if (c > 0) {
                cur = remove_pos_swap(cur, r, threads);
                frontier -= 1;
            }
        }
    }

    while (frontier > 0) {
        cur = remove_last(cur, threads);
        frontier -= 1;
    }

    Poly result{};
    for (size_t i = 0; i < cur.states.size(); ++i) {
        if (cur.states[i].size == 0) {
            result = cur.polys[i];
            break;
        }
    }
    return result;
}

static int64_t pow_mod(int64_t a, int64_t e) {
    int64_t res = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1) res = (res * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return res;
}

static int64_t eval_poly(const Poly& coeffs, int max_deg, int64_t n) {
    int64_t res = 0;
    int64_t x = n % MOD;
    for (int d = max_deg; d >= 0; --d) {
        res = (res * x + coeffs[d]) % MOD;
    }
    return res;
}

static int64_t lagrange_eval(const vector<int64_t>& y,
                             int64_t x,
                             const vector<int64_t>& fact,
                             const vector<int64_t>& inv_fact) {
    int m = static_cast<int>(y.size()) - 1;
    if (x <= m) return y[static_cast<size_t>(x)];

    vector<int64_t> prefix(static_cast<size_t>(m) + 1);
    vector<int64_t> suffix(static_cast<size_t>(m) + 1);
    prefix[0] = 1;
    for (int i = 0; i < m; ++i) {
        int64_t term = x - i;
        if (term < 0) term += MOD;
        prefix[static_cast<size_t>(i + 1)] = (prefix[static_cast<size_t>(i)] * term) % MOD;
    }
    suffix[static_cast<size_t>(m)] = 1;
    for (int i = m - 1; i >= 0; --i) {
        int64_t term = x - (i + 1);
        if (term < 0) term += MOD;
        suffix[static_cast<size_t>(i)] = (suffix[static_cast<size_t>(i + 1)] * term) % MOD;
    }

    int64_t res = 0;
    for (int i = 0; i <= m; ++i) {
        int64_t num = prefix[static_cast<size_t>(i)] * suffix[static_cast<size_t>(i)] % MOD;
        int64_t denom = inv_fact[static_cast<size_t>(i)] * inv_fact[static_cast<size_t>(m - i)] % MOD;
        int64_t term = y[static_cast<size_t>(i)] * num % MOD * denom % MOD;
        if ((m - i) & 1) term = MOD - term;
        res += term;
        if (res >= MOD) res -= MOD;
    }
    return res;
}

static int64_t sum_pows(int d, int64_t n, const vector<int64_t>& fact, const vector<int64_t>& inv_fact) {
    if (n <= 0) return 0;
    if (d == 0) return n % MOD;
    int m = d + 1;
    vector<int64_t> y(static_cast<size_t>(m) + 1);
    y[0] = 0;
    for (int i = 1; i <= m; ++i) {
        y[static_cast<size_t>(i)] = (y[static_cast<size_t>(i - 1)] + pow_mod(i, d)) % MOD;
    }
    int64_t x = n % MOD;
    return lagrange_eval(y, x, fact, inv_fact);
}

static int64_t compute_S(const Poly& coeffs, int max_deg, int64_t n,
                         const vector<int64_t>& fact, const vector<int64_t>& inv_fact) {
    int64_t res = 0;
    for (int d = 0; d <= max_deg; ++d) {
        if (coeffs[d] == 0) continue;
        int64_t sum = sum_pows(d, n, fact, inv_fact);
        res = (res + coeffs[d] * sum) % MOD;
    }
    return res;
}

static bool run_validation(int threads, const vector<int64_t>& fact, const vector<int64_t>& inv_fact) {
    bool ok = true;

    {
        Poly p = chromatic_poly_grid(2, 2, threads);
        int64_t f3 = eval_poly(p, 4, 3);
        int64_t f20 = eval_poly(p, 4, 20);
        if (f3 != 18) {
            cerr << "Validation failed: F(2,2,3) got " << f3 << ", expected 18\n";
            ok = false;
        }
        if (f20 != 130340) {
            cerr << "Validation failed: F(2,2,20) got " << f20 << ", expected 130340\n";
            ok = false;
        }
    }

    {
        Poly p = chromatic_poly_grid(3, 4, threads);
        int64_t f6 = eval_poly(p, 12, 6);
        if (f6 != 102923670) {
            cerr << "Validation failed: F(3,4,6) got " << f6 << ", expected 102923670\n";
            ok = false;
        }
    }

    {
        Poly p = chromatic_poly_grid(4, 4, threads);
        int64_t s15 = compute_S(p, 16, 15, fact, inv_fact);
        if (s15 != 325951319) {
            cerr << "Validation failed: S(4,4,15) got " << s15 << ", expected 325951319\n";
            ok = false;
        }
    }

    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    bool validate = true;
    unsigned hw = thread::hardware_concurrency();
    int threads = hw ? static_cast<int>(hw) : 1;
    threads = max(1, min(threads, 8));

    for (int i = 1; i < argc; ++i) {
        string arg(argv[i]);
        if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = max(1, stoi(argv[++i]));
        }
    }

    vector<int64_t> fact(static_cast<size_t>(MAX_DEG) + 3, 1);
    vector<int64_t> inv_fact(static_cast<size_t>(MAX_DEG) + 3, 1);
    for (size_t i = 1; i < fact.size(); ++i) fact[i] = fact[i - 1] * static_cast<int64_t>(i) % MOD;
    inv_fact.back() = pow_mod(fact.back(), MOD - 2);
    for (size_t i = fact.size() - 1; i > 0; --i) {
        inv_fact[i - 1] = inv_fact[i] * static_cast<int64_t>(i) % MOD;
    }

    if (validate && !run_validation(min(threads, 2), fact, inv_fact)) {
        return 1;
    }

    constexpr int R = 9;
    constexpr int C = 10;
    constexpr int64_t N = 1112131415LL;
    Poly poly = chromatic_poly_grid(R, C, threads);
    int64_t answer = compute_S(poly, R * C, N, fact, inv_fact);
    cout << answer << "\n";
    return 0;
}
