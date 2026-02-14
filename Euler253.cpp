#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

struct VecHash {
    std::size_t operator()(const std::vector<int>& v) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (int x : v) {
            h ^= static_cast<std::size_t>(x + 0x9e3779b97f4a7c15ULL);
            h *= 1099511628211ULL;
        }
        return h;
    }
};

using State = std::vector<int>;

struct Solver {
    std::unordered_map<State, std::vector<long double>, VecHash> memo;
    std::unordered_map<State, std::vector<std::pair<State, int>>, VecHash> trans;

    static State make_next(const State& s, int remove_len, int add1, int add2 = 0) {
        State t = s;
        auto it = std::lower_bound(t.begin(), t.end(), remove_len);
        t.erase(it);
        if (add1 > 0) {
            t.insert(std::lower_bound(t.begin(), t.end(), add1), add1);
        }
        if (add2 > 0) {
            t.insert(std::lower_bound(t.begin(), t.end(), add2), add2);
        }
        return t;
    }

    const std::vector<std::pair<State, int>>& get_transitions(const State& s) {
        auto it = trans.find(s);
        if (it != trans.end()) return it->second;

        std::unordered_map<State, int, VecHash> acc;
        const int blocks = static_cast<int>(s.size());
        int i = 0;
        while (i < blocks) {
            int j = i;
            while (j < blocks && s[j] == s[i]) ++j;
            const int len = s[i];
            const int cnt = j - i;

            if (len == 1) {
                State t = make_next(s, 1, 0, 0);
                acc[t] += cnt;
            } else {
                State endpoint = make_next(s, len, len - 1, 0);
                acc[endpoint] += 2 * cnt;
                for (int a = 1; a <= len - 2; ++a) {
                    int b = len - 1 - a;
                    State split = make_next(s, len, a, b);
                    acc[split] += cnt;
                }
            }
            i = j;
        }

        std::vector<std::pair<State, int>> out;
        out.reserve(acc.size());
        for (auto& kv : acc) out.push_back(std::move(kv));
        auto inserted = trans.emplace(s, std::move(out));
        return inserted.first->second;
    }

    long double expected_max(const State& s, int current_max) {
        if (s.empty()) return static_cast<long double>(current_max);
        auto& row = memo[s];
        if (row.empty()) {
            row.assign(41, std::numeric_limits<long double>::quiet_NaN());
        }
        if (row[current_max] == row[current_max]) return row[current_max];

        int remaining = 0;
        for (int x : s) remaining += x;
        const auto& tr = get_transitions(s);

        long double ans = 0;
        for (const auto& [ns, w] : tr) {
            int next_max = std::max(current_max, static_cast<int>(ns.size()));
            ans += (static_cast<long double>(w) / remaining) * expected_max(ns, next_max);
        }
        row[current_max] = ans;
        return ans;
    }
};

static int max_segments_of_perm(const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    std::vector<char> on(n + 2, 0);
    int seg = 0;
    int best = 0;
    for (int pos : p) {
        on[pos] = 1;
        const bool l = on[pos - 1];
        const bool r = on[pos + 1];
        if (l && r) --seg;
        else if (!l && !r) ++seg;
        if (seg > best) best = seg;
    }
    return best;
}

static long double brute_expectation(int n) {
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 1);
    long double total = 0;
    long double cnt = 0;
    do {
        total += max_segments_of_perm(p);
        cnt += 1;
    } while (std::next_permutation(p.begin(), p.end()));
    return total / cnt;
}

int main() {
    {
        Solver t1;
        long double dp6 = t1.expected_max(State{6}, 1);
        long double br6 = brute_expectation(6);
        assert(std::fabsl(dp6 - br6) < 1e-12L);

        Solver t2;
        long double dp7 = t2.expected_max(State{7}, 1);
        long double br7 = brute_expectation(7);
        assert(std::fabsl(dp7 - br7) < 1e-12L);
    }

    Solver solver;
    long double ans = solver.expected_max(State{40}, 1);
    std::cout << std::fixed << std::setprecision(6) << static_cast<double>(ans) << '\n';
    return 0;
}
