#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1'000'000'007LL;
constexpr int kMainFact = 200;
constexpr u64 kMainN = 1'000'000'000'000ULL;
constexpr int kSmallLimit = 1'000'000;

constexpr int kCheckFact1 = 3;
constexpr u64 kCheckN1 = 100ULL;
constexpr i64 kCheckExpected1 = 3398LL;
constexpr int kCheckFact2 = 4;
constexpr u64 kCheckN2 = 1'000'000ULL;
constexpr i64 kCheckExpected2 = 268'882'292LL;

inline i64 mod_norm(i64 x) {
    x %= kMod;
    if (x < 0) x += kMod;
    return x;
}

inline i64 mod_mul(const i64 a, const i64 b) {
    return static_cast<i64>((static_cast<i128>(a) * b) % kMod);
}

i64 mod_pow(i64 base, i64 exp) {
    i64 result = 1 % kMod;
    base = mod_norm(base);
    while (exp > 0) {
        if (exp & 1LL) result = mod_mul(result, base);
        base = mod_mul(base, base);
        exp >>= 1LL;
    }
    return result;
}

std::vector<int> make_primes(const int n) {
    std::vector<unsigned char> mark(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    for (int p = 2; p <= n; ++p) {
        if (mark[static_cast<std::size_t>(p)] != 0) continue;
        primes.push_back(p);
        for (int q = p; q <= n; q += p) {
            mark[static_cast<std::size_t>(q)] = 1;
        }
    }
    return primes;
}

int exponent_in_factorial(int p, int n) {
    int e = 0;
    while (n > 0) {
        n /= p;
        e += n;
    }
    return e;
}

struct U64Hash {
    std::size_t operator()(u64 x) const noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
        return static_cast<std::size_t>(x ^ (x >> 31U));
    }
};

std::vector<i64> small_prefix_tau;
std::unordered_map<u64, i64, U64Hash> large_prefix_tau_cache;

void build_small_prefix_tau() {
    small_prefix_tau.assign(static_cast<std::size_t>(kSmallLimit), 0);
    for (int d = 1; d < kSmallLimit; ++d) {
        for (int m = d; m < kSmallLimit; m += d) {
            ++small_prefix_tau[static_cast<std::size_t>(m)];
        }
    }
    for (int i = 1; i < kSmallLimit; ++i) {
        i64 v = small_prefix_tau[static_cast<std::size_t>(i)] +
                small_prefix_tau[static_cast<std::size_t>(i - 1)];
        if (v >= kMod) v -= kMod;
        small_prefix_tau[static_cast<std::size_t>(i)] = v;
    }
}

i64 prefix_tau(const u64 n) {
    if (n < static_cast<u64>(kSmallLimit)) {
        return small_prefix_tau[static_cast<std::size_t>(n)];
    }
    const auto it = large_prefix_tau_cache.find(n);
    if (it != large_prefix_tau_cache.end()) {
        return it->second;
    }

    i128 sum = 0;
    u64 x = 1;
    while (x * x <= n) {
        sum += static_cast<i128>(n / x);
        ++x;
    }
    const i128 val = 2 * sum - static_cast<i128>(x - 1) * static_cast<i128>(x - 1);
    const i64 out = mod_norm(static_cast<i64>(val % kMod));
    large_prefix_tau_cache.emplace(n, out);
    return out;
}

struct Solver608 {
    std::vector<int> primes;
    std::vector<i64> ds;
    u64 N = 0;
    i64 ans = 0;

    static i64 tri(const i64 x) { return x * (x + 1LL) / 2LL; }

    void dfs(const u64 acc, const int start_idx, const i64 d, const int sign) {
        const i64 term = mod_mul(d, prefix_tau(N / acc));
        ans = (sign > 0) ? mod_norm(ans + term) : mod_norm(ans - term);

        for (int j = start_idx; j < static_cast<int>(primes.size()); ++j) {
            const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(j)]);
            if (acc > N / p) break;
            dfs(acc * p, j + 1, mod_mul(d, ds[static_cast<std::size_t>(j)]), -sign);
        }
    }

    i64 compute(const int fact_n, const u64 n) {
        N = n;
        ans = 0;
        primes = make_primes(fact_n);
        ds.assign(primes.size(), 0);

        i64 k = 1;
        for (std::size_t i = 0; i < primes.size(); ++i) {
            const int e = exponent_in_factorial(primes[i], fact_n);
            const i64 t_e = tri(e) % kMod;
            const i64 t_ep1 = tri(e + 1) % kMod;
            ds[i] = mod_mul(t_e, mod_pow(t_ep1, kMod - 2));
            k = mod_mul(k, t_ep1);
        }

        dfs(1ULL, 0, k, +1);
        return ans;
    }
};

}  // namespace

int main() {
    build_small_prefix_tau();
    large_prefix_tau_cache.reserve(1 << 16);

    Solver608 solver;

    assert(solver.compute(kCheckFact1, kCheckN1) == mod_norm(kCheckExpected1));
    assert(solver.compute(kCheckFact2, kCheckN2) == mod_norm(kCheckExpected2));

    std::cout << solver.compute(kMainFact, kMainN) << '\n';
    return 0;
}
