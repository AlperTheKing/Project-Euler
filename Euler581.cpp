#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 581: 47-smooth triangular numbers
//
// T(n) = n(n+1)/2 is 47-smooth  <=>  its prime factors are all <= 47.
// Since gcd(n, n+1)=1, after dividing by 2 the two coprime factors must each be 47-smooth:
// - If n is even:  T(n) = (n/2) * (n+1)  =>  n/2 and n+1 are 47-smooth.
// - If n is odd :  T(n) = n * ((n+1)/2) =>  n and (n+1)/2 are 47-smooth.
//
// Let x be the 47-smooth factor:
// - n = 2x      requires 2x+1 to be 47-smooth.
// - n = 2x-1    requires 2x-1 to be 47-smooth.
//
// So we enumerate all 47-smooth x up to a safe bound and test 2x±1.

using u64 = std::uint64_t;

static constexpr u64 LIMIT_X = 10'000'000'000'000ULL;  // 1e13

static const int PRIMES[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

static inline bool is_47_smooth(u64 x) {
    for (int p : PRIMES) {
        while (x % (u64)p == 0) x /= (u64)p;
    }
    return x == 1;
}

static u64 sum_n = 0;

static void dfs(int idx, u64 cur) {
    if (idx == (int)(sizeof(PRIMES) / sizeof(PRIMES[0]))) {
        const u64 x = cur;
        const u64 y_plus = 2 * x + 1;
        const u64 y_minus = 2 * x - 1;

        if (is_47_smooth(y_plus)) sum_n += 2 * x;
        if (is_47_smooth(y_minus)) sum_n += 2 * x - 1;
        return;
    }
    const u64 p = (u64)PRIMES[idx];
    u64 v = cur;
    while (v <= LIMIT_X) {
        dfs(idx + 1, v);
        if (v > LIMIT_X / p) break;
        v *= p;
    }
}

int main() {
    dfs(0, 1);
    std::cout << sum_n << "\n";
    return 0;
}
