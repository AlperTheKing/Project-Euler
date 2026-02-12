#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <iostream>

using i64 = std::int64_t;
using i128 = __int128_t;

struct Interval {
    int l;
    int u;
};

struct Option {
    int bin;
    int l;
    int u;
};

struct VecHash {
    std::size_t operator()(const std::vector<int>& v) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (int x : v) {
            h ^= static_cast<std::size_t>(x + 0x9e3779b9);
            h *= 1099511628211ULL;
        }
        return h;
    }
};

static bool interval_less(const Interval& a, const Interval& b) {
    if (a.l != b.l) {
        return a.l < b.l;
    }
    return a.u < b.u;
}

static int lcm_upto(int n) {
    i64 d = 1;
    for (int i = 1; i <= n; ++i) {
        d = std::lcm(d, static_cast<i64>(i));
    }
    return static_cast<int>(d);
}

class Solver {
public:
    explicit Solver(int target) : target_(target), d_(lcm_upto(target)) {}

    int denom() const {
        return d_;
    }

    i64 solve_min() {
        std::vector<Interval> init{{0, d_}};
        return solve_min_state(init);
    }

    bool can_choose_all() {
        std::vector<Interval> init{{0, d_}};
        return feasible_state(init);
    }

private:
    int target_;
    int d_;
    std::unordered_map<std::vector<int>, i64, VecHash> memo_min_;
    std::unordered_map<std::vector<int>, bool, VecHash> memo_feasible_;

    std::vector<int> pack(const std::vector<Interval>& intervals) const {
        std::vector<int> key;
        key.reserve(intervals.size() * 2);
        for (const auto& in : intervals) {
            key.push_back(in.l);
            key.push_back(in.u);
        }
        return key;
    }

    i64 solve_min_state(const std::vector<Interval>& intervals) {
        const std::vector<int> key = pack(intervals);
        auto it = memo_min_.find(key);
        if (it != memo_min_.end()) {
            return it->second;
        }

        const int n = static_cast<int>(intervals.size());
        if (n == target_) {
            i64 sum = 0;
            for (const auto& in : intervals) {
                sum += in.l;
            }
            memo_min_.emplace(key, sum);
            return sum;
        }

        const int m = n + 1;
        std::vector<std::pair<int, int>> bins(m);
        for (int k = 0; k < m; ++k) {
            bins[k] = {k * d_ / m, (k + 1) * d_ / m};
        }

        std::vector<std::vector<Option>> options(n);
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < m; ++k) {
                const int nl = std::max(intervals[i].l, bins[k].first);
                const int nu = std::min(intervals[i].u, bins[k].second);
                if (nl < nu) {
                    options[i].push_back({k, nl, nu});
                }
            }
            if (options[i].empty()) {
                static constexpr i64 INF = (1LL << 60);
                memo_min_.emplace(key, INF);
                return INF;
            }
        }

        std::vector<int> order(n);
        for (int i = 0; i < n; ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return options[a].size() < options[b].size();
        });

        std::vector<char> used(m, 0);
        std::vector<Interval> assigned(n);
        i64 best = (1LL << 60);

        std::function<void(int)> dfs = [&](int pos) {
            if (pos == n) {
                int hole = -1;
                for (int k = 0; k < m; ++k) {
                    if (!used[k]) {
                        hole = k;
                        break;
                    }
                }
                if (hole < 0) {
                    return;
                }
                std::vector<Interval> child = assigned;
                child.push_back({bins[hole].first, bins[hole].second});
                std::sort(child.begin(), child.end(), interval_less);

                const i64 cand = solve_min_state(child);
                if (cand < best) {
                    best = cand;
                }
                return;
            }

            const int i = order[pos];
            for (const auto& op : options[i]) {
                if (used[op.bin]) {
                    continue;
                }
                used[op.bin] = 1;
                assigned[i] = {op.l, op.u};
                dfs(pos + 1);
                used[op.bin] = 0;
            }
        };

        dfs(0);
        memo_min_.emplace(key, best);
        return best;
    }

    bool feasible_state(const std::vector<Interval>& intervals) {
        const std::vector<int> key = pack(intervals);
        auto it = memo_feasible_.find(key);
        if (it != memo_feasible_.end()) {
            return it->second;
        }

        const int n = static_cast<int>(intervals.size());
        if (n == target_) {
            memo_feasible_.emplace(key, true);
            return true;
        }

        const int m = n + 1;
        std::vector<std::pair<int, int>> bins(m);
        for (int k = 0; k < m; ++k) {
            bins[k] = {k * d_ / m, (k + 1) * d_ / m};
        }

        std::vector<std::vector<Option>> options(n);
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < m; ++k) {
                const int nl = std::max(intervals[i].l, bins[k].first);
                const int nu = std::min(intervals[i].u, bins[k].second);
                if (nl < nu) {
                    options[i].push_back({k, nl, nu});
                }
            }
            if (options[i].empty()) {
                memo_feasible_.emplace(key, false);
                return false;
            }
        }

        std::vector<int> order(n);
        for (int i = 0; i < n; ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return options[a].size() < options[b].size();
        });

        std::vector<char> used(m, 0);
        std::vector<Interval> assigned(n);
        bool ok = false;

        std::function<void(int)> dfs = [&](int pos) {
            if (ok) {
                return;
            }
            if (pos == n) {
                int hole = -1;
                for (int k = 0; k < m; ++k) {
                    if (!used[k]) {
                        hole = k;
                        break;
                    }
                }
                if (hole < 0) {
                    return;
                }
                std::vector<Interval> child = assigned;
                child.push_back({bins[hole].first, bins[hole].second});
                std::sort(child.begin(), child.end(), interval_less);
                if (feasible_state(child)) {
                    ok = true;
                }
                return;
            }

            const int i = order[pos];
            for (const auto& op : options[i]) {
                if (used[op.bin]) {
                    continue;
                }
                used[op.bin] = 1;
                assigned[i] = {op.l, op.u};
                dfs(pos + 1);
                used[op.bin] = 0;
                if (ok) {
                    return;
                }
            }
        };

        dfs(0);
        memo_feasible_.emplace(key, ok);
        return ok;
    }
};

static i64 brute_min(int target) {
    const int d = lcm_upto(target);

    std::function<i64(const std::vector<Interval>&)> dfs = [&](const std::vector<Interval>& intervals) -> i64 {
        const int n = static_cast<int>(intervals.size());
        if (n == target) {
            i64 sum = 0;
            for (const auto& in : intervals) {
                sum += in.l;
            }
            return sum;
        }

        const int m = n + 1;
        std::vector<std::pair<int, int>> bins(m);
        for (int k = 0; k < m; ++k) {
            bins[k] = {k * d / m, (k + 1) * d / m};
        }

        std::vector<std::vector<Option>> options(n);
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < m; ++k) {
                const int nl = std::max(intervals[i].l, bins[k].first);
                const int nu = std::min(intervals[i].u, bins[k].second);
                if (nl < nu) {
                    options[i].push_back({k, nl, nu});
                }
            }
            if (options[i].empty()) {
                return (1LL << 60);
            }
        }

        std::vector<char> used(m, 0);
        std::vector<Interval> assigned(n);
        i64 best = (1LL << 60);

        std::function<void(int)> assign = [&](int idx) {
            if (idx == n) {
                int hole = -1;
                for (int k = 0; k < m; ++k) {
                    if (!used[k]) {
                        hole = k;
                        break;
                    }
                }
                if (hole < 0) {
                    return;
                }
                std::vector<Interval> child = assigned;
                child.push_back({bins[hole].first, bins[hole].second});
                std::sort(child.begin(), child.end(), interval_less);
                best = std::min(best, dfs(child));
                return;
            }

            for (const auto& op : options[idx]) {
                if (used[op.bin]) {
                    continue;
                }
                used[op.bin] = 1;
                assigned[idx] = {op.l, op.u};
                assign(idx + 1);
                used[op.bin] = 0;
            }
        };

        assign(0);
        return best;
    };

    std::vector<Interval> init{{0, d}};
    return dfs(init);
}

static std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static std::string format_rounded_12(i64 numerator, i64 denominator) {
    static constexpr i128 SCALE = 1'000'000'000'000LL;
    const i128 num = static_cast<i128>(numerator);
    const i128 den = static_cast<i128>(denominator);

    const i128 rounded = (num * SCALE * 2 + den) / (den * 2);
    const i128 integer_part = rounded / SCALE;
    const i128 frac_part = rounded % SCALE;

    std::string frac = to_string_i128(frac_part);
    if (frac.size() < 12) {
        frac = std::string(12 - frac.size(), '0') + frac;
    }

    return to_string_i128(integer_part) + "." + frac;
}

int main() {
    {
        Solver check4(4);
        const i64 f4 = check4.solve_min();
        assert(2 * f4 == 3LL * check4.denom());
    }

    {
        Solver fast6(6);
        const i64 dp6 = fast6.solve_min();
        const i64 brute6 = brute_min(6);
        assert(dp6 == brute6);
    }

    {
        Solver feasible17(17);
        assert(feasible17.can_choose_all());
        Solver feasible18(18);
        assert(!feasible18.can_choose_all());
    }

    Solver solver17(17);
    const i64 num = solver17.solve_min();
    std::cout << format_rounded_12(num, solver17.denom()) << '\n';
    return 0;
}
