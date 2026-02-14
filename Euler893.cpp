#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;

static constexpr int DIGIT_COST[10] = {6, 2, 5, 5, 4, 5, 6, 3, 7, 6};

static bool improve_mul(std::vector<int>& m) {
    const int N = static_cast<int>(m.size()) - 1;
    bool changed = false;
    for (int x = 1; x <= N; ++x) {
        const int cx = m[static_cast<std::size_t>(x)];
        const int maxy = N / x;
        for (int y = x; y <= maxy; ++y) {
            const int prod = x * y;
            const int cand = cx + m[static_cast<std::size_t>(y)] + 2;
            if (cand < m[static_cast<std::size_t>(prod)]) {
                m[static_cast<std::size_t>(prod)] = cand;
                changed = true;
            }
        }
    }
    return changed;
}

static bool improve_add(std::vector<int>& m) {
    const int N = static_cast<int>(m.size()) - 1;
    const int maxm = *std::max_element(m.begin() + 1, m.end());
    std::vector<std::vector<int>> by_cost(static_cast<std::size_t>(maxm + 1));
    for (int n = 1; n <= N; ++n) by_cost[static_cast<std::size_t>(m[static_cast<std::size_t>(n)])].push_back(n);

    bool changed = false;
    for (int ca = 0; ca <= maxm; ++ca) {
        const auto& A = by_cost[static_cast<std::size_t>(ca)];
        if (A.empty()) continue;
        for (int cb = ca; cb <= maxm; ++cb) {
            const auto& B = by_cost[static_cast<std::size_t>(cb)];
            if (B.empty()) continue;
            const int cand = ca + cb + 2;
            if (cand >= maxm) break;

            if (ca == cb) {
                const int as = static_cast<int>(A.size());
                for (int i = 0; i < as; ++i) {
                    const int x = A[static_cast<std::size_t>(i)];
                    const int limit = N - x;
                    auto end_it = std::upper_bound(A.begin() + i, A.end(), limit);
                    for (auto it = A.begin() + i; it != end_it; ++it) {
                        const int s = x + *it;
                        if (cand < m[static_cast<std::size_t>(s)]) {
                            m[static_cast<std::size_t>(s)] = cand;
                            changed = true;
                        }
                    }
                }
            } else {
                const std::vector<int>* P = &A;
                const std::vector<int>* Q = &B;
                if (P->size() > Q->size()) std::swap(P, Q);
                for (int x : *P) {
                    const int limit = N - x;
                    auto end_it = std::upper_bound(Q->begin(), Q->end(), limit);
                    for (auto it = Q->begin(); it != end_it; ++it) {
                        const int s = x + *it;
                        if (cand < m[static_cast<std::size_t>(s)]) {
                            m[static_cast<std::size_t>(s)] = cand;
                            changed = true;
                        }
                    }
                }
            }
        }
    }
    return changed;
}

static u64 solve_fast(int N) {
    std::vector<int> m(static_cast<std::size_t>(N + 1), 0);
    for (int n = 1; n <= N; ++n) {
        m[static_cast<std::size_t>(n)] = m[static_cast<std::size_t>(n / 10)] + DIGIT_COST[n % 10];
    }
    while (improve_mul(m)) {}
    while (improve_add(m)) {}
    u64 sum = 0;
    for (int n = 1; n <= N; ++n) sum += static_cast<u64>(m[static_cast<std::size_t>(n)]);
    return sum;
}

static u64 solve_exact_small(int N) {
    std::vector<int> A(static_cast<std::size_t>(N + 1), 0);
    for (int n = 1; n <= N; ++n) {
        A[static_cast<std::size_t>(n)] = A[static_cast<std::size_t>(n / 10)] + DIGIT_COST[n % 10];
    }
    for (int n = 2; n <= N; ++n) {
        int best = A[static_cast<std::size_t>(n)];
        int root = static_cast<int>(std::sqrt(static_cast<double>(n)));
        for (int d = 2; d <= root; ++d) {
            if (n % d == 0) {
                int q = n / d;
                int cand = A[static_cast<std::size_t>(d)] + A[static_cast<std::size_t>(q)] + 2;
                if (cand < best) best = cand;
            }
        }
        A[static_cast<std::size_t>(n)] = best;
    }

    std::vector<int> M = A;
    for (int n = 1; n <= N; ++n) {
        int best = M[static_cast<std::size_t>(n)];
        for (int a = 1; a < n; ++a) {
            int cand = M[static_cast<std::size_t>(a)] + A[static_cast<std::size_t>(n - a)] + 2;
            if (cand < best) best = cand;
        }
        M[static_cast<std::size_t>(n)] = best;
    }
    u64 sum = 0;
    for (int n = 1; n <= N; ++n) sum += static_cast<u64>(M[static_cast<std::size_t>(n)]);
    return sum;
}

int main() {
    assert(solve_fast(100) == 916ULL);
    assert(solve_fast(2000) == solve_exact_small(2000));
    assert(solve_fast(5000) == solve_exact_small(5000));
    std::cout << solve_fast(1'000'000) << '\n';
    return 0;
}
