#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

using u64 = std::uint64_t;

static constexpr int DIGIT_COST[10] = {6, 2, 5, 5, 4, 5, 6, 3, 7, 6};

static std::vector<int> compute_A(int N) {
    std::vector<int> A(static_cast<std::size_t>(N + 1), 0);
    for (int n = 1; n <= N; ++n) {
        A[static_cast<std::size_t>(n)] = A[static_cast<std::size_t>(n / 10)] + DIGIT_COST[n % 10];
    }

    for (int n = 2; n <= N; ++n) {
        int root = static_cast<int>(std::sqrt(static_cast<double>(n)));
        for (int d = 2; d <= root; ++d) {
            if (n % d != 0) continue;
            int q = n / d;
            int cand = A[static_cast<std::size_t>(d)] + A[static_cast<std::size_t>(q)] + 2;
            if (cand < A[static_cast<std::size_t>(n)]) A[static_cast<std::size_t>(n)] = cand;
        }
    }

    return A;
}

static u64 T_fast(int N) {
    std::vector<int> A = compute_A(N);
    std::vector<int> M = A;

    int maxA = 0;
    for (int n = 1; n <= N; ++n) maxA = std::max(maxA, A[static_cast<std::size_t>(n)]);

    std::vector<std::vector<int>> by_cost(static_cast<std::size_t>(maxA + 1));
    for (int n = 1; n <= N; ++n) by_cost[static_cast<std::size_t>(A[static_cast<std::size_t>(n)])].push_back(n);

    int maxM = maxA;
    for (int c = 0; c <= maxA; ++c) {
        const auto& v = by_cost[static_cast<std::size_t>(c)];
        for (int a : v) {
            int ca = c + 2;
            int s = a + 1;
            int maxb = N - a;
            for (int b = 1; b <= maxb; ++b, ++s) {
                int cand = ca + A[static_cast<std::size_t>(b)];
                if (cand < M[static_cast<std::size_t>(s)]) M[static_cast<std::size_t>(s)] = cand;
            }
        }

        maxM = 0;
        for (int n = 1; n <= N; ++n) maxM = std::max(maxM, M[static_cast<std::size_t>(n)]);
        if (maxM <= 2 * c + 3) break;
    }

    u64 sum = 0;
    for (int n = 1; n <= N; ++n) sum += static_cast<u64>(M[static_cast<std::size_t>(n)]);
    return sum;
}

static u64 T_exact_small(int N) {
    std::vector<int> A = compute_A(N);
    std::vector<int> M = A;

    for (int n = 1; n <= N; ++n) {
        int best = A[static_cast<std::size_t>(n)];
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
    assert(T_fast(100) == 916ULL);
    assert(T_fast(2000) == T_exact_small(2000));
    assert(T_fast(5000) == T_exact_small(5000));

    std::cout << T_fast(1'000'000) << '\n';
    return 0;
}
