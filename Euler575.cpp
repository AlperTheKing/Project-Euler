#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

// Project Euler 575: Wandering Robot
//
// Two possible transition rules on an N x N grid:
// (i) Each available option (stay + each adjacent move) has equal probability:
//     P(stay)=1/(d+1), P(move to each neighbor)=1/(d+1).
// (ii) Fixed 50% probability to stay, other 50% split among neighbors:
//     P(stay)=1/2, P(move to each neighbor)=1/(2d).
//
// Both are random walks on an undirected weighted graph, hence stationary probability
// is proportional to the node's weighted degree:
// (i) weight = d+1  (self-loop weight 1, each edge weight 1)
// (ii) weight = d    (self-loop weight d, each edge weight 1 -> total 2d, proportional to d)
//
// The final answer is the average of the two stationary probabilities (coin flip).

using u64 = std::uint64_t;

static long double probability_square_room(int N) {
    const u64 NN = (u64)N * (u64)N;

    const u64 corners = 4;
    const u64 edges = 4ULL * (N - 2);
    const u64 interior = (u64)(N - 2) * (u64)(N - 2);

    const u64 total_deg = 2 * corners + 3 * edges + 4 * interior;
    const u64 total_deg_plus_1 = total_deg + NN;  // sum(d+1) = sum(d) + number of nodes

    u64 sum_deg_sq = 0;
    for (int k = 1; k <= N; ++k) {
        const u64 x = (u64)k * (u64)k;  // room number (1-indexed)
        const u64 r = (x - 1) / (u64)N + 1;  // 1..N
        const u64 c = (x - 1) % (u64)N + 1;  // 1..N

        const bool top = (r == 1), bottom = (r == (u64)N);
        const bool left = (c == 1), right = (c == (u64)N);
        const bool is_corner = (top || bottom) && (left || right);
        const bool is_edge = (top || bottom || left || right);

        const u64 d = is_corner ? 2ULL : (is_edge ? 3ULL : 4ULL);
        sum_deg_sq += d;
    }

    const u64 count_sq = (u64)N;
    const u64 sum_deg_plus_1_sq = sum_deg_sq + count_sq;

    const long double p_i = (long double)sum_deg_plus_1_sq / (long double)total_deg_plus_1;
    const long double p_ii = (long double)sum_deg_sq / (long double)total_deg;
    return 0.5L * (p_i + p_ii);
}

int main() {
    {
        const long double p5 = probability_square_room(5);
        const long double expected = 0.177976190476L;
        if (fabsl(p5 - expected) > 5e-13L) {  // tolerate 0.5 ulp at 12 d.p.
            std::cerr << "Validation failed for 5x5 grid\n";
            return 1;
        }
    }

    const long double ans = probability_square_room(1000);
    std::cout << std::fixed << std::setprecision(12) << (double)ans << "\n";
    return 0;
}

