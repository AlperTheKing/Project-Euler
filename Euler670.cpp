#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1'000'004'321ULL;
constexpr int COLORS = 4;

struct State {
    int rt;
    int ct;
    int rb;
    int cb;
    int prev_vert;

    bool operator==(const State& other) const {
        return rt == other.rt && ct == other.ct && rb == other.rb && cb == other.cb &&
               prev_vert == other.prev_vert;
    }

    bool operator<(const State& other) const {
        if (rt != other.rt) return rt < other.rt;
        if (ct != other.ct) return ct < other.ct;
        if (rb != other.rb) return rb < other.rb;
        if (cb != other.cb) return cb < other.cb;
        return prev_vert < other.prev_vert;
    }
};

struct StateHash {
    std::size_t operator()(const State& s) const {
        std::size_t h = static_cast<std::size_t>(s.rt);
        h = h * 1315423911u + static_cast<std::size_t>(s.ct + 7);
        h = h * 1315423911u + static_cast<std::size_t>(s.rb + 11);
        h = h * 1315423911u + static_cast<std::size_t>(s.cb + 13);
        h = h * 1315423911u + static_cast<std::size_t>(s.prev_vert + 17);
        return h;
    }
};

std::vector<State> initial_states() {
    std::vector<State> init;
    init.reserve(4 + 3 * 3 * 4 * 3);

    for (int v = 0; v < COLORS; ++v) {
        init.push_back({0, v, 0, v, 1});
    }

    for (int lt = 1; lt <= 3; ++lt) {
        for (int lb = 1; lb <= 3; ++lb) {
            for (int a = 0; a < COLORS; ++a) {
                for (int b = 0; b < COLORS; ++b) {
                    if (a == b) {
                        continue;
                    }
                    init.push_back({lt - 1, a, lb - 1, b, 0});
                }
            }
        }
    }

    return init;
}

std::vector<State> next_states(const State& s) {
    std::vector<State> out;

    const bool top_change = (s.rt == 0);
    const bool bottom_change = (s.rb == 0);

    if (s.rt == 0 && s.rb == 0) {
        for (int v = 0; v < COLORS; ++v) {
            if (v == s.ct || v == s.cb) {
                continue;
            }
            out.push_back({0, v, 0, v, 1});
        }
    }

    std::vector<std::pair<int, int>> top_options;
    std::vector<std::pair<int, int>> bottom_options;

    if (s.rt > 0) {
        top_options.push_back({s.rt - 1, s.ct});
    } else {
        for (int len = 1; len <= 3; ++len) {
            for (int color = 0; color < COLORS; ++color) {
                if (color == s.ct) {
                    continue;
                }
                top_options.push_back({len - 1, color});
            }
        }
    }

    if (s.rb > 0) {
        bottom_options.push_back({s.rb - 1, s.cb});
    } else {
        for (int len = 1; len <= 3; ++len) {
            for (int color = 0; color < COLORS; ++color) {
                if (color == s.cb) {
                    continue;
                }
                bottom_options.push_back({len - 1, color});
            }
        }
    }

    for (const auto& top : top_options) {
        for (const auto& bottom : bottom_options) {
            const int nt_rt = top.first;
            const int nt_ct = top.second;
            const int nt_rb = bottom.first;
            const int nt_cb = bottom.second;

            if (nt_ct == nt_cb) {
                continue;
            }

            if (top_change && bottom_change && s.prev_vert == 0) {
                continue;
            }

            out.push_back({nt_rt, nt_ct, nt_rb, nt_cb, 0});
        }
    }

    return out;
}

struct Matrix {
    int n;
    std::vector<u64> a;

    explicit Matrix(int n_) : n(n_), a(static_cast<std::size_t>(n_) * static_cast<std::size_t>(n_), 0ULL) {}

    u64& at(int i, int j) { return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)]; }
    u64 at(int i, int j) const {
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
    }
};

Matrix multiply(const Matrix& x, const Matrix& y) {
    const int n = x.n;
    Matrix z(n);

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            const u64 xik = x.at(i, k);
            if (xik == 0ULL) {
                continue;
            }
            for (int j = 0; j < n; ++j) {
                const u64 ykj = y.at(k, j);
                if (ykj == 0ULL) {
                    continue;
                }
                u64& cell = z.at(i, j);
                cell = static_cast<u64>((static_cast<u128>(cell) + static_cast<u128>(xik) * static_cast<u128>(ykj)) % MOD);
            }
        }
    }

    return z;
}

std::vector<u64> multiply(const std::vector<u64>& v, const Matrix& m) {
    const int n = m.n;
    std::vector<u64> out(static_cast<std::size_t>(n), 0ULL);

    for (int i = 0; i < n; ++i) {
        const u64 vi = v[static_cast<std::size_t>(i)];
        if (vi == 0ULL) {
            continue;
        }
        for (int j = 0; j < n; ++j) {
            const u64 mij = m.at(i, j);
            if (mij == 0ULL) {
                continue;
            }
            u64& cell = out[static_cast<std::size_t>(j)];
            cell = static_cast<u64>((static_cast<u128>(cell) + static_cast<u128>(vi) * static_cast<u128>(mij)) % MOD);
        }
    }

    return out;
}

class Solver {
public:
    Solver() {
        init_states_ = initial_states();
        build_reachable();
        build_transition();
        build_initial_vector();
    }

    u64 count(u64 n) const {
        assert(n >= 1ULL);

        std::vector<u64> vec = init_vector_;
        Matrix p = transition_;
        u64 exp = n - 1ULL;

        while (exp > 0ULL) {
            if ((exp & 1ULL) != 0ULL) {
                vec = multiply(vec, p);
            }
            exp >>= 1ULL;
            if (exp > 0ULL) {
                p = multiply(p, p);
            }
        }

        u64 ans = 0ULL;
        for (int i = 0; i < static_cast<int>(states_.size()); ++i) {
            const State& s = states_[static_cast<std::size_t>(i)];
            if (s.rt == 0 && s.rb == 0) {
                ans += vec[static_cast<std::size_t>(i)];
                if (ans >= MOD) {
                    ans -= MOD;
                }
            }
        }
        return ans;
    }

private:
    std::vector<State> init_states_;
    std::vector<State> states_;
    std::unordered_map<State, int, StateHash> index_;
    Matrix transition_{1};
    std::vector<u64> init_vector_;

    void build_reachable() {
        std::unordered_map<State, int, StateHash> seen;
        std::queue<State> q;

        for (const State& s : init_states_) {
            if (seen.emplace(s, 1).second) {
                q.push(s);
            }
        }

        while (!q.empty()) {
            const State s = q.front();
            q.pop();

            const std::vector<State> ns = next_states(s);
            for (const State& t : ns) {
                if (seen.emplace(t, 1).second) {
                    q.push(t);
                }
            }
        }

        states_.reserve(seen.size());
        for (const auto& kv : seen) {
            states_.push_back(kv.first);
        }
        std::sort(states_.begin(), states_.end());

        index_.reserve(states_.size() * 2);
        for (int i = 0; i < static_cast<int>(states_.size()); ++i) {
            index_[states_[static_cast<std::size_t>(i)]] = i;
        }
    }

    void build_transition() {
        const int n = static_cast<int>(states_.size());
        transition_ = Matrix(n);

        for (int i = 0; i < n; ++i) {
            const std::vector<State> ns = next_states(states_[static_cast<std::size_t>(i)]);
            for (const State& t : ns) {
                const int j = index_.at(t);
                u64& cell = transition_.at(i, j);
                ++cell;
                if (cell >= MOD) {
                    cell -= MOD;
                }
            }
        }
    }

    void build_initial_vector() {
        init_vector_.assign(states_.size(), 0ULL);
        for (const State& s : init_states_) {
            const int i = index_.at(s);
            u64& cell = init_vector_[static_cast<std::size_t>(i)];
            ++cell;
            if (cell >= MOD) {
                cell -= MOD;
            }
        }
    }
};

}  // namespace

int main() {
    Solver solver;

    assert(solver.count(2ULL) == 120ULL);
    assert(solver.count(5ULL) == 45'876ULL);
    assert(solver.count(100ULL) == 53'275'818ULL);

    std::cout << solver.count(10'000'000'000'000'000ULL) << "\n";
    return 0;
}
