#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

struct Outcome {
    int x, y;
    long double p;
};

static long double solve_optimal_expected(int cards, const std::vector<Outcome> &outs) {
    const int S = 1 << cards;
    const int goal = S - 1;

    struct Act {
        uint16_t idx;
    };
    std::vector<std::array<Act, 3>> pol((size_t)S * outs.size());
    auto pol_at = [&](int s, int oi) -> std::array<Act, 3> & { return pol[(size_t)s * outs.size() + oi]; };

    std::vector<long double> V(S, 0);
    V[goal] = 0;

    for (int s = 0; s < S; ++s) {
        for (int oi = 0; oi < (int)outs.size(); ++oi) {
            const int x = outs[oi].x, y = outs[oi].y;
            const int a1 = x - 1;
            const int a2 = y - 1;
            const int a3 = x + y - 1;
            pol_at(s, oi) = {Act{(uint16_t)a1}, Act{(uint16_t)a2}, Act{(uint16_t)a3}};
        }
    }

    auto eval_policy = [&](int iters, long double eps) {
        for (int it = 0; it < iters; ++it) {
            long double delta = 0;
            for (int s = 0; s < S; ++s) {
                if (s == goal) continue;
                long double ex = 0;
                for (int oi = 0; oi < (int)outs.size(); ++oi) {
                    const auto &acts = pol_at(s, oi);
                    const int nxt = s ^ (1 << acts[0].idx);
                    ex += outs[oi].p * V[nxt];
                }
                const long double nv = 1 + ex;
                delta = std::max(delta, fabsl(nv - V[s]));
                V[s] = nv;
            }
            if (delta < eps) break;
        }
    };

    auto improve_policy = [&]() -> bool {
        bool changed = false;
        for (int s = 0; s < S; ++s) {
            if (s == goal) continue;
            for (int oi = 0; oi < (int)outs.size(); ++oi) {
                const auto &acts = pol_at(s, oi);
                int best_idx = 0;
                long double best_v = 1e300L;
                for (int k = 0; k < 3; ++k) {
                    const int nxt = s ^ (1 << acts[k].idx);
                    const long double v = V[nxt];
                    if (v < best_v) {
                        best_v = v;
                        best_idx = k;
                    }
                }
                if (best_idx != 0) {
                    auto &w = pol_at(s, oi);
                    std::swap(w[0], w[best_idx]);
                    changed = true;
                }
            }
        }
        return changed;
    };

    for (int round = 0; round < 50; ++round) {
        eval_policy(10'000, 1e-15L);
        if (!improve_policy()) break;
    }

    eval_policy(50'000, 1e-18L);
    return V[0];
}

int main() {
    {
        std::vector<Outcome> coin;
        for (int x : {1, 2}) {
            for (int y : {1, 2}) coin.push_back({x, y, 0.25L});
        }
        const long double e = solve_optimal_expected(4, coin);
        assert(fabsl(e - 5.673651L) < 1e-6L);
    }

    std::vector<Outcome> dice;
    dice.reserve(36);
    for (int x = 1; x <= 6; ++x) {
        for (int y = 1; y <= 6; ++y) dice.push_back({x, y, 1.0L / 36.0L});
    }
    const long double e = solve_optimal_expected(12, dice);

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6) << (double)e << "\n";
    return 0;
}
