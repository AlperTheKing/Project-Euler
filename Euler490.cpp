#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr uint64_t MOD = 1000000000ULL;

struct State {
    int k = 0;
    std::array<int, 4> deg{0, 0, 0, 0};
    std::array<int, 4> comp{0, 0, 0, 0};
};

struct StateHash {
    size_t operator()(const State& s) const noexcept {
        size_t h = static_cast<size_t>(s.k);
        for (int i = 0; i < s.k; ++i) {
            h = h * 3 + static_cast<size_t>(s.deg[i]);
        }
        for (int i = 0; i < s.k; ++i) {
            h = h * 3 + static_cast<size_t>(s.comp[i]);
        }
        return h;
    }
};

struct StateEq {
    bool operator()(const State& a, const State& b) const noexcept {
        if (a.k != b.k) return false;
        for (int i = 0; i < a.k; ++i) {
            if (a.deg[i] != b.deg[i]) return false;
        }
        for (int i = 0; i < a.k; ++i) {
            if (a.comp[i] != b.comp[i]) return false;
        }
        return true;
    }
};

State canonicalize(State s) {
    int map[5];
    for (int& v : map) v = -1;
    int next = 0;
    for (int i = 0; i < s.k; ++i) {
        int c = s.comp[i];
        if (map[c] == -1) map[c] = next++;
        s.comp[i] = map[c];
    }
    return s;
}

bool comps_all_same(const State& s) {
    for (int i = 1; i < s.k; ++i) {
        if (s.comp[i] != s.comp[0]) return false;
    }
    return true;
}

std::vector<std::pair<State, uint64_t>> transitions(const State& s, bool remove, bool remove_is_one) {
    std::vector<std::pair<State, uint64_t>> out;
    int k = s.k;
    std::vector<std::vector<int>> subsets;
    subsets.push_back({});
    for (int i = 0; i < k; ++i) subsets.push_back({i});
    for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) {
            subsets.push_back({i, j});
        }
    }

    auto add_state = [&](const State& ns) {
        for (auto& entry : out) {
            if (StateEq{}(entry.first, ns)) {
                entry.second += 1;
                return;
            }
        }
        out.push_back({ns, 1});
    };

    for (const auto& subset : subsets) {
        if (subset.size() == 2 && s.comp[subset[0]] == s.comp[subset[1]]) continue;

        State ns;
        ns.k = k + 1;
        ns.deg = s.deg;
        ns.comp = s.comp;

        bool ok = true;
        for (int idx : subset) {
            if (ns.deg[idx] >= 2) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;

        for (int idx : subset) ns.deg[idx] += 1;
        int new_deg = static_cast<int>(subset.size());
        int new_comp = 0;

        if (subset.empty()) {
            int maxc = -1;
            for (int i = 0; i < k; ++i) maxc = std::max(maxc, ns.comp[i]);
            new_comp = maxc + 1;
        } else if (subset.size() == 1) {
            new_comp = ns.comp[subset[0]];
        } else {
            int comp_a = ns.comp[subset[0]];
            int comp_b = ns.comp[subset[1]];
            new_comp = std::min(comp_a, comp_b);
            for (int i = 0; i < k; ++i) {
                if (ns.comp[i] == comp_a || ns.comp[i] == comp_b) ns.comp[i] = new_comp;
            }
        }

        ns.deg[k] = new_deg;
        ns.comp[k] = new_comp;

        if (remove) {
            int leaving_deg = ns.deg[0];
            if (remove_is_one) {
                if (leaving_deg != 1) continue;
            } else {
                if (leaving_deg != 2) continue;
            }
            int leaving_comp = ns.comp[0];
            for (int i = 1; i <= k; ++i) {
                ns.deg[i - 1] = ns.deg[i];
                ns.comp[i - 1] = ns.comp[i];
            }
            ns.k = k;
            bool present = false;
            for (int i = 0; i < ns.k; ++i) {
                if (ns.comp[i] == leaving_comp) {
                    present = true;
                    break;
                }
            }
            if (!present) continue;
        }

        ns = canonicalize(ns);
        add_state(ns);
    }

    return out;
}

bool is_final_small(const State& s, int n) {
    if (!comps_all_same(s)) return false;
    if (n == 1) {
        return s.k == 1 && s.deg[0] == 0;
    }
    if (n == 2) {
        return s.k == 2 && s.deg[0] == 1 && s.deg[1] == 1;
    }
    if (n == 3) {
        return s.k == 3 && s.deg[0] == 1 && s.deg[1] == 2 && s.deg[2] == 1;
    }
    return s.k == 3 && s.deg[0] == 2 && s.deg[1] == 2 && s.deg[2] == 1;
}

bool is_final_generic(const State& s) {
    return s.k == 3 && comps_all_same(s) && s.deg[0] == 2 && s.deg[1] == 2 && s.deg[2] == 1;
}

std::vector<uint64_t> compute_small_f(int n_max) {
    std::unordered_map<State, uint64_t, StateHash, StateEq> dp;
    State start;
    start.k = 1;
    start.deg = {0, 0, 0, 0};
    start.comp = {0, 0, 0, 0};
    dp[start] = 1;

    std::vector<uint64_t> f(n_max + 1, 0);
    for (int t = 1; t <= n_max; ++t) {
        uint64_t total = 0;
        for (const auto& kv : dp) {
            if (is_final_small(kv.first, t)) total += kv.second;
        }
        f[t] = total;
        if (t == n_max) break;
        bool remove = (t + 1 >= 4);
        bool remove_is_one = remove && (t - 2 == 1);
        std::unordered_map<State, uint64_t, StateHash, StateEq> next;
        for (const auto& kv : dp) {
            for (const auto& nxt : transitions(kv.first, remove, remove_is_one)) {
                next[nxt.first] += kv.second * nxt.second;
            }
        }
        dp.swap(next);
    }
    return f;
}

struct Model {
    int d = 0;
    std::vector<State> states;
    std::unordered_map<State, int, StateHash, StateEq> index;
    std::vector<std::vector<uint64_t>> A;
    std::vector<uint64_t> v4;
    std::vector<uint64_t> c;
};

Model build_model() {
    std::unordered_map<State, uint64_t, StateHash, StateEq> dp;
    State start;
    start.k = 1;
    start.deg = {0, 0, 0, 0};
    start.comp = {0, 0, 0, 0};
    dp[start] = 1;

    for (int t = 1; t <= 3; ++t) {
        bool remove = (t + 1 >= 4);
        bool remove_is_one = remove && (t - 2 == 1);
        std::unordered_map<State, uint64_t, StateHash, StateEq> next;
        for (const auto& kv : dp) {
            for (const auto& nxt : transitions(kv.first, remove, remove_is_one)) {
                uint64_t add = (kv.second * nxt.second) % MOD;
                uint64_t& slot = next[nxt.first];
                slot += add;
                if (slot >= MOD) slot %= MOD;
            }
        }
        dp.swap(next);
    }

    Model model;
    std::queue<State> q;
    for (const auto& kv : dp) {
        if (model.index.emplace(kv.first, model.d).second) {
            model.states.push_back(kv.first);
            q.push(kv.first);
            ++model.d;
        }
    }

    while (!q.empty()) {
        State cur = q.front();
        q.pop();
        for (const auto& nxt : transitions(cur, true, false)) {
            if (model.index.emplace(nxt.first, model.d).second) {
                model.states.push_back(nxt.first);
                q.push(nxt.first);
                ++model.d;
            }
        }
    }

    model.A.assign(model.d, std::vector<uint64_t>(model.d, 0));
    for (int j = 0; j < model.d; ++j) {
        const State& cur = model.states[j];
        for (const auto& nxt : transitions(cur, true, false)) {
            int i = model.index[nxt.first];
            model.A[i][j] = (model.A[i][j] + nxt.second) % MOD;
        }
    }

    model.v4.assign(model.d, 0);
    for (const auto& kv : dp) {
        int id = model.index[kv.first];
        model.v4[id] = (model.v4[id] + kv.second) % MOD;
    }

    model.c.assign(model.d, 0);
    for (int i = 0; i < model.d; ++i) {
        if (is_final_generic(model.states[i])) model.c[i] = 1;
    }

    return model;
}

using Matrix = std::vector<std::vector<uint64_t>>;

Matrix mul_mat(const Matrix& A, const Matrix& B) {
    int d = static_cast<int>(A.size());
    Matrix C(d, std::vector<uint64_t>(d, 0));
    for (int i = 0; i < d; ++i) {
        for (int k = 0; k < d; ++k) {
            if (A[i][k] == 0) continue;
            for (int j = 0; j < d; ++j) {
                __int128 term = static_cast<__int128>(A[i][k]) * B[k][j];
                C[i][j] = static_cast<uint64_t>((C[i][j] + term) % MOD);
            }
        }
    }
    return C;
}

std::vector<uint64_t> mul_mat_vec(const Matrix& A, const std::vector<uint64_t>& v) {
    int d = static_cast<int>(A.size());
    std::vector<uint64_t> out(d, 0);
    for (int i = 0; i < d; ++i) {
        __int128 sum = 0;
        for (int j = 0; j < d; ++j) {
            sum += static_cast<__int128>(A[i][j]) * v[j];
        }
        out[i] = static_cast<uint64_t>(sum % MOD);
    }
    return out;
}

int idx3(int i, int j, int k, int d) {
    return (i * d + j) * d + k;
}

std::vector<uint64_t> transform_tensor(const std::vector<uint64_t>& T, const Matrix& A) {
    int d = static_cast<int>(A.size());
    int size = d * d * d;
    std::vector<uint64_t> U(size, 0), V(size, 0), W(size, 0);

    for (int a = 0; a < d; ++a) {
        for (int j = 0; j < d; ++j) {
            for (int k = 0; k < d; ++k) {
                __int128 sum = 0;
                for (int i = 0; i < d; ++i) {
                    sum += static_cast<__int128>(A[i][a]) * T[idx3(i, j, k, d)];
                }
                U[idx3(a, j, k, d)] = static_cast<uint64_t>(sum % MOD);
            }
        }
    }

    for (int a = 0; a < d; ++a) {
        for (int b = 0; b < d; ++b) {
            for (int k = 0; k < d; ++k) {
                __int128 sum = 0;
                for (int j = 0; j < d; ++j) {
                    sum += static_cast<__int128>(A[j][b]) * U[idx3(a, j, k, d)];
                }
                V[idx3(a, b, k, d)] = static_cast<uint64_t>(sum % MOD);
            }
        }
    }

    for (int a = 0; a < d; ++a) {
        for (int b = 0; b < d; ++b) {
            for (int c = 0; c < d; ++c) {
                __int128 sum = 0;
                for (int k = 0; k < d; ++k) {
                    sum += static_cast<__int128>(A[k][c]) * V[idx3(a, b, k, d)];
                }
                W[idx3(a, b, c, d)] = static_cast<uint64_t>(sum % MOD);
            }
        }
    }

    return W;
}

uint64_t eval_tensor(const std::vector<uint64_t>& T, const std::vector<uint64_t>& v) {
    int d = static_cast<int>(v.size());
    uint64_t total = 0;
    for (int i = 0; i < d; ++i) {
        if (v[i] == 0) continue;
        for (int j = 0; j < d; ++j) {
            if (v[j] == 0) continue;
            uint64_t vij = (v[i] * v[j]) % MOD;
            for (int k = 0; k < d; ++k) {
                if (v[k] == 0) continue;
                uint64_t t = T[idx3(i, j, k, d)];
                if (t == 0) continue;
                __int128 term = static_cast<__int128>(t) * vij % MOD;
                term = term * v[k] % MOD;
                total += static_cast<uint64_t>(term);
                if (total >= MOD) total %= MOD;
            }
        }
    }
    return total % MOD;
}

struct Powers {
    std::vector<Matrix> A_pow;
    std::vector<std::vector<uint64_t>> T_pow;
};

Powers build_powers(const Matrix& A, const std::vector<uint64_t>& c, int max_bits) {
    int d = static_cast<int>(A.size());
    Powers powers;
    powers.A_pow.resize(max_bits);
    powers.T_pow.resize(max_bits);

    powers.A_pow[0] = A;
    std::vector<uint64_t> T1(d * d * d, 0);
    for (int i = 0; i < d; ++i) {
        if (c[i] == 0) continue;
        for (int j = 0; j < d; ++j) {
            if (c[j] == 0) continue;
            uint64_t cij = (c[i] * c[j]) % MOD;
            for (int k = 0; k < d; ++k) {
                if (c[k] == 0) continue;
                T1[idx3(i, j, k, d)] = (cij * c[k]) % MOD;
            }
        }
    }
    powers.T_pow[0] = std::move(T1);

    for (int bit = 1; bit < max_bits; ++bit) {
        powers.A_pow[bit] = mul_mat(powers.A_pow[bit - 1], powers.A_pow[bit - 1]);
        std::vector<uint64_t> transformed = transform_tensor(powers.T_pow[bit - 1], powers.A_pow[bit - 1]);
        std::vector<uint64_t> merged(transformed.size(), 0);
        for (size_t i = 0; i < transformed.size(); ++i) {
            uint64_t val = powers.T_pow[bit - 1][i] + transformed[i];
            if (val >= MOD) val %= MOD;
            merged[i] = val;
        }
        powers.T_pow[bit] = std::move(merged);
    }

    return powers;
}

uint64_t sum_from_4(uint64_t L, const Model& model, const Powers& powers) {
    if (L < 4) return 0;
    uint64_t len = L - 3;
    std::vector<uint64_t> v = model.v4;
    uint64_t sum = 0;
    int bit = 0;
    while (len > 0) {
        if (len & 1ULL) {
            sum = (sum + eval_tensor(powers.T_pow[bit], v)) % MOD;
            v = mul_mat_vec(powers.A_pow[bit], v);
        }
        len >>= 1ULL;
        ++bit;
    }
    return sum;
}

uint64_t compute_S_mod(uint64_t L, const Model& model, const Powers& powers) {
    uint64_t total = 0;
    if (L >= 1) total = (total + 1) % MOD;
    if (L >= 2) total = (total + 1) % MOD;
    if (L >= 3) total = (total + 1) % MOD;
    if (L >= 4) total = (total + sum_from_4(L, model, powers)) % MOD;
    return total;
}

bool validate(const Model& model, const Powers& powers) {
    auto f = compute_small_f(40);
    if (f[6] != 14) {
        std::cerr << "Validation failed: f(6) = " << f[6] << "\n";
        return false;
    }
    if (f[10] != 254) {
        std::cerr << "Validation failed: f(10) = " << f[10] << "\n";
        return false;
    }
    if (f[40] != 1439682432976ULL) {
        std::cerr << "Validation failed: f(40) = " << f[40] << "\n";
        return false;
    }

    __int128 sum10 = 0;
    for (int i = 1; i <= 10; ++i) {
        __int128 v = f[i];
        sum10 += v * v * v;
    }
    if (static_cast<uint64_t>(sum10) != 18230635ULL) {
        std::cerr << "Validation failed: S(10) = " << static_cast<uint64_t>(sum10) << "\n";
        return false;
    }

    __int128 sum20 = 0;
    for (int i = 1; i <= 20; ++i) {
        __int128 v = f[i];
        sum20 += v * v * v;
    }
    if (static_cast<uint64_t>(sum20) != 104207881192114219ULL) {
        std::cerr << "Validation failed: S(20) = " << static_cast<uint64_t>(sum20) << "\n";
        return false;
    }

    uint64_t s1000 = compute_S_mod(1000, model, powers);
    if (s1000 != 225031475ULL) {
        std::cerr << "Validation failed: S(1000) mod 1e9 = " << s1000 << "\n";
        return false;
    }

    uint64_t s1e6 = compute_S_mod(1000000ULL, model, powers);
    if (s1e6 != 363486179ULL) {
        std::cerr << "Validation failed: S(1e6) mod 1e9 = " << s1e6 << "\n";
        return false;
    }

    return true;
}

} // namespace

int main(int argc, char** argv) {
    uint64_t L = 100000000000000ULL;
    if (argc > 1) {
        L = std::stoull(argv[1]);
    }

    Model model = build_model();
    uint64_t max_L = L < 1000000ULL ? 1000000ULL : L;
    uint64_t max_len = (max_L >= 4) ? (max_L - 3) : 1;
    int max_bits = 0;
    while (max_bits < 63 && (1ULL << max_bits) <= max_len) ++max_bits;
    if (max_bits == 0) max_bits = 1;
    Powers powers = build_powers(model.A, model.c, max_bits);

    if (!validate(model, powers)) return 1;

    uint64_t answer = compute_S_mod(L, model, powers);
    std::cout << answer << '\n';
    return 0;
}
