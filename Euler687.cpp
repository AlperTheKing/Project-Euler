#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <functional>

namespace {

using ld = long double;

struct ProbDP {
    std::unordered_map<std::uint32_t, ld> memo;

    static std::uint32_t encode(int f, int c1, int c2, int c3, int c4, int prev) {
        const int p = (prev < 0) ? 0 : prev;
        return static_cast<std::uint32_t>(f) |
               (static_cast<std::uint32_t>(c1) << 6U) |
               (static_cast<std::uint32_t>(c2) << 10U) |
               (static_cast<std::uint32_t>(c3) << 14U) |
               (static_cast<std::uint32_t>(c4) << 18U) |
               (static_cast<std::uint32_t>(p) << 22U);
    }

    ld dfs(int f, int c1, int c2, int c3, int c4, int prev) {
        const int rem = f + c1 + 2 * c2 + 3 * c3 + 4 * c4;
        if (rem == 0) {
            return 1.0L;
        }

        const std::uint32_t key = encode(f, c1, c2, c3, c4, prev);
        const auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        ld out = 0.0L;

        if (f > 0) {
            out += (static_cast<ld>(f) / static_cast<ld>(rem)) *
                   dfs(f - 1, c1, c2, c3, c4, -1);
        }

        auto add_rank = [&](int r, int count, int nc1, int nc2, int nc3, int nc4) {
            if (count == 0) {
                return;
            }
            int avail = count;
            if (prev == r) {
                --avail;
            }
            if (avail <= 0) {
                return;
            }
            const int np = (r - 1 == 0) ? -1 : (r - 1);
            const ld prob = static_cast<ld>(avail * r) / static_cast<ld>(rem);
            out += prob * dfs(f, nc1, nc2, nc3, nc4, np);
        };

        add_rank(1, c1, c1 - 1, c2, c3, c4);
        add_rank(2, c2, c1 + 1, c2 - 1, c3, c4);
        add_rank(3, c3, c1, c2 + 1, c3 - 1, c4);
        add_rank(4, c4, c1, c2, c3 + 1, c4 - 1);

        memo.emplace(key, out);
        return out;
    }

    ld prob_all_perfect_for_subset_size(int j) {
        memo.clear();
        memo.reserve(1 << 18);
        return dfs(52 - 4 * j, 0, 0, 0, j, -1);
    }
};

}  // namespace

int main() {
    std::array<std::array<int, 14>, 14> C{};
    for (int n = 0; n <= 13; ++n) {
        C[static_cast<std::size_t>(n)][0] = 1;
        C[static_cast<std::size_t>(n)][static_cast<std::size_t>(n)] = 1;
        for (int k = 1; k < n; ++k) {
            C[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)] =
                C[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(k - 1)] +
                C[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(k)];
        }
    }

    ProbDP solver;

    std::array<ld, 14> p{};
    p[0] = 1.0L;
    for (int j = 1; j <= 13; ++j) {
        p[static_cast<std::size_t>(j)] = solver.prob_all_perfect_for_subset_size(j);
    }

    std::array<ld, 14> M{};
    for (int j = 0; j <= 13; ++j) {
        M[static_cast<std::size_t>(j)] =
            static_cast<ld>(C[13][j]) * p[static_cast<std::size_t>(j)];
    }

    std::array<ld, 14> P{};
    for (int k = 13; k >= 0; --k) {
        ld v = M[static_cast<std::size_t>(k)];
        for (int t = k + 1; t <= 13; ++t) {
            v -= static_cast<ld>(C[static_cast<std::size_t>(t)][static_cast<std::size_t>(k)]) *
                 P[static_cast<std::size_t>(t)];
        }
        P[static_cast<std::size_t>(k)] = v;
    }

    ld expect = 0.0L;
    for (int k = 0; k <= 13; ++k) {
        expect += static_cast<ld>(k) * P[static_cast<std::size_t>(k)];
    }
    assert(std::fabsl(expect - (4324.0L / 425.0L)) < 1e-12L);

    const std::array<int, 6> primes{{2, 3, 5, 7, 11, 13}};
    ld ans = 0.0L;
    for (int pval : primes) {
        ans += P[static_cast<std::size_t>(pval)];
    }

    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(ans) << "\n";
    return 0;
}
