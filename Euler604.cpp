#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

// Project Euler 604: Maximum lattice points on a strictly convex increasing graph in an N x N square.
// After reducing to primitive step-vectors (a,b) with a,b>0 and gcd(a,b)=1, each step costs (a+b) in L1.
// Taking all steps with a+b=s gives exactly phi(s) choices and total L1 cost s*phi(s), so the L1-greedy
// optimum is found from prefix sums of s*phi(s), with a single parity/balancing exception when L1 is tight.

using u64 = std::uint64_t;

static bool has_balanced_triple_even2(u64 s) {
    // s == 2*n with n odd. Need 3 totatives a,b,c of s summing to 3s/2.
    const u64 n = s / 2;
    if (n < 7) return false; // covers (6,10) style small failures quickly

    const u64 LIM = 2000; // small even offsets are enough for n in our range
    for (u64 d1 = 2; d1 <= LIM && d1 < n; d1 += 2) {
        if (std::gcd(d1, n) != 1) continue;
        for (u64 d2 = d1 + 2; d2 <= LIM && d1 + d2 < n; d2 += 2) {
            if (std::gcd(d2, n) != 1) continue;
            if (std::gcd(d1 + d2, n) != 1) continue;
            return true;
        }
    }
    return false;
}

static u64 F(u64 N, const std::vector<u64>& pref_phi, const std::vector<u64>& pref_cost) {
    const u64 target = 2 * N; // total L1 budget

    const auto it = std::upper_bound(pref_cost.begin(), pref_cost.end(), target);
    const u64 m = (it == pref_cost.begin()) ? 0 : (u64)(it - pref_cost.begin() - 1);

    const u64 base_cnt = (m >= 2) ? pref_phi[m] : 0;   // sum_{s=2..m} phi(s)
    const u64 base_cost = (m >= 2) ? pref_cost[m] : 0; // sum_{s=2..m} s*phi(s)

    const u64 slack = target - base_cost;
    const u64 s = m + 1; // next L1 layer
    const u64 t = (s > 0) ? (slack / s) : 0;
    const u64 r = (s > 0) ? (slack % s) : 0;

    u64 segments = base_cnt + t;

    // If r==0 and t is odd, L1 is saturated and we must have exact x=y=N; this requires a balanced odd subset
    // in the last layer. For s divisible by 4 it's impossible by parity; for s==2*n (n odd) it's equivalent to
    // existence of a balanced triple (then pad with complementary pairs).
    if (r == 0 && (t & 1) && t > 0) {
        bool ok = false;
        if (t == 1) {
            ok = false;
        } else if ((s % 4) == 0) {
            ok = false;
        } else if ((s % 4) == 2) {
            ok = has_balanced_triple_even2(s);
        }
        if (!ok) --segments;
    }

    return segments + 1; // lattice points = segments + start point
}

int main() {
    const u64 N_max = 1000000000000000000ULL;
    const u64 target_max = 2 * N_max;

    // Sieve totients until prefix L1 cost covers 2*N_max.
    u64 limit = 3000000;
    std::vector<int> phi;
    std::vector<u64> pref_phi;
    std::vector<u64> pref_cost;
    for (;;) {
        phi.assign((size_t)limit + 1, 0);
        for (u64 i = 0; i <= limit; ++i) phi[i] = (int)i;
        for (u64 p = 2; p <= limit; ++p) {
            if ((u64)phi[p] != p) continue;
            for (u64 j = p; j <= limit; j += p) phi[j] -= phi[j] / (int)p;
        }

        pref_phi.assign((size_t)limit + 1, 0);
        pref_cost.assign((size_t)limit + 1, 0);
        for (u64 i = 2; i <= limit; ++i) {
            pref_phi[i] = pref_phi[i - 1] + (u64)phi[i];
            pref_cost[i] = pref_cost[i - 1] + i * (u64)phi[i];
        }

        if (pref_cost[limit] >= target_max) break;
        limit = (u64)(limit * 1.4) + 1000;
    }

    // Statement validations.
    if (F(1, pref_phi, pref_cost) != 2) {
        std::cerr << "Validation failed: F(1)\n";
        return 1;
    }
    if (F(3, pref_phi, pref_cost) != 3) {
        std::cerr << "Validation failed: F(3)\n";
        return 1;
    }
    if (F(9, pref_phi, pref_cost) != 6) {
        std::cerr << "Validation failed: F(9)\n";
        return 1;
    }
    if (F(11, pref_phi, pref_cost) != 7) {
        std::cerr << "Validation failed: F(11)\n";
        return 1;
    }
    if (F(100, pref_phi, pref_cost) != 30) {
        std::cerr << "Validation failed: F(100)\n";
        return 1;
    }
    if (F(50000, pref_phi, pref_cost) != 1898) {
        std::cerr << "Validation failed: F(50000)\n";
        return 1;
    }

    std::cout << F(N_max, pref_phi, pref_cost) << "\n";
    return 0;
}

