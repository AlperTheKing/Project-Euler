#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using boost::multiprecision::uint256_t;

constexpr u64 kMod = 1'000'000'007ULL;

struct Factor {
    int p;
    int a;
};

struct Elem {
    uint256_t v;
    u64 w;
};

uint256_t pow10_u256(int e) {
    uint256_t x = 1;
    for (int i = 0; i < e; ++i) x *= 10;
    return x;
}

std::vector<int> primes_upto(int n) {
    std::vector<char> is_prime(static_cast<std::size_t>(n + 1), 1);
    if (n >= 0) is_prime[0] = 0;
    if (n >= 1) is_prime[1] = 0;
    for (int p = 2; (long long)p * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) continue;
        for (int m = p * p; m <= n; m += p) is_prime[static_cast<std::size_t>(m)] = 0;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

std::vector<Factor> factorial_factors(int n) {
    std::vector<Factor> out;
    const auto primes = primes_upto(n);
    for (int p : primes) {
        int a = 0;
        for (int x = n; x; x /= p) a += x / p;
        out.push_back(Factor{p, a});
    }
    return out;
}

u64 mod_mul(u64 a, u64 b) { return (a * b) % kMod; }

void gen_divs(const std::vector<Factor>& fac, const int idx, const uint256_t& cap, uint256_t v, u64 w,
              std::vector<Elem>& out) {
    if (idx == (int)fac.size()) {
        out.push_back(Elem{v, w});
        return;
    }
    const int p = fac[(std::size_t)idx].p;
    const int a = fac[(std::size_t)idx].a;
    const u64 negp = (kMod - static_cast<u64>(p) % kMod) % kMod;
    uint256_t vv = v;
    u64 ww = w;
    for (int e = 0; e <= a; ++e) {
        gen_divs(fac, idx + 1, cap, vv, ww, out);
        if (e == a) break;
        if (vv > cap / static_cast<unsigned>(p)) break;
        vv *= static_cast<unsigned>(p);
        ww = mod_mul(ww, negp);
    }
}

u64 prefix_sum(const uint256_t& X, const std::vector<Elem>& A, const std::vector<Elem>& B,
               const std::vector<u64>& prefB) {
    std::size_t j = B.size();
    u64 ans = 0;
    for (const auto& a : A) {
        if (a.v > X) break;
        const uint256_t lim = X / a.v;
        while (j > 0 && B[j - 1].v > lim) --j;
        ans += mod_mul(a.w, prefB[j]);
        ans %= kMod;
    }
    return ans;
}

u64 S_factorial_bounded(const int n, const uint256_t& L, const uint256_t& H) {
    const std::vector<Factor> fac = factorial_factors(n);
    const int split = std::min<int>(5, static_cast<int>(fac.size()));
    const std::vector<Factor> fa(fac.begin(), fac.begin() + split);
    const std::vector<Factor> fb(fac.begin() + split, fac.end());

    u64 cntA = 1;
    for (const auto& f : fa) cntA *= static_cast<u64>(f.a + 1);
    u64 cntB = 1;
    for (const auto& f : fb) cntB *= static_cast<u64>(f.a + 1);

    std::vector<Elem> A;
    std::vector<Elem> B;
    A.reserve(static_cast<std::size_t>(cntA));
    B.reserve(static_cast<std::size_t>(cntB));

    gen_divs(fa, 0, H, 1, 1, A);
    gen_divs(fb, 0, H, 1, 1, B);

    std::sort(A.begin(), A.end(), [](const Elem& x, const Elem& y) { return x.v < y.v; });
    std::sort(B.begin(), B.end(), [](const Elem& x, const Elem& y) { return x.v < y.v; });

    std::vector<u64> prefB(B.size() + 1, 0);
    for (std::size_t i = 0; i < B.size(); ++i) prefB[i + 1] = (prefB[i] + B[i].w) % kMod;

    const u64 sumH = prefix_sum(H, A, B, prefB);
    const uint256_t Lm1 = L - 1;
    const u64 sumL = (L == 0 || L == 1) ? 0 : prefix_sum(Lm1, A, B, prefB);
    return (sumH + kMod - sumL) % kMod;
}

}  // namespace

int main() {
    assert(S_factorial_bounded(10, 100, 1000) == 1457);
    assert(S_factorial_bounded(15, 1000, 100000) == (kMod - 107'974ULL));
    assert(S_factorial_bounded(30, 100'000'000, 1'000'000'000'000ULL) == (9766732243224ULL % kMod));

    const uint256_t L = pow10_u256(20);
    const uint256_t H = pow10_u256(60);
    std::cout << S_factorial_bounded(70, L, H) << "\n";
    return 0;
}
