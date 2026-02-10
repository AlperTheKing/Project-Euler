#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 619: vectors of squarefree kernels over GF(2); count is 2^(n-rank)-1.

using u64 = std::uint64_t;

static std::vector<int> sieve_spf(int n) {
    std::vector<int> spf(n + 1, 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            if ((long long)i * i <= n) {
                for (long long j = 1LL * i * i; j <= n; j += i)
                    if (spf[(int)j] == 0) spf[(int)j] = i;
            }
        }
    }
    spf[1] = 1;
    return spf;
}

static std::vector<int> list_primes_from_spf(const std::vector<int> &spf) {
    std::vector<int> primes;
    for (int i = 2; i < (int)spf.size(); ++i)
        if (spf[i] == i) primes.push_back(i);
    return primes;
}

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e) {
        if (e & 1ULL) r = (u64)((__int128)r * a % mod);
        a = (u64)((__int128)a * a % mod);
        e >>= 1ULL;
    }
    return r;
}

static void xor_inplace(std::vector<int> &v, const std::vector<int> &b, std::vector<int> &tmp) {
    tmp.clear();
    tmp.reserve(v.size() + b.size());
    std::size_t i = 0, j = 0;
    while (i < v.size() || j < b.size()) {
        if (j == b.size() || (i < v.size() && v[i] < b[j])) {
            tmp.push_back(v[i++]);
        } else if (i == v.size() || b[j] < v[i]) {
            tmp.push_back(b[j++]);
        } else {
            ++i;
            ++j;
        }
    }
    v.swap(tmp);
}

static int rank_squarefree_vectors(int a, int b, const std::vector<int> &spf, const std::vector<int> &prime_id,
                                   int num_primes) {
    std::vector<std::vector<int>> basis((std::size_t)num_primes);
    std::vector<int> tmp;
    int rank = 0;
    for (int x = a; x <= b; ++x) {
        int n = x;
        std::vector<int> v;
        while (n > 1) {
            int p = spf[n];
            int cnt = 0;
            while (n % p == 0) {
                n /= p;
                cnt ^= 1;
            }
            if (cnt) v.push_back(prime_id[p]);
        }
        while (!v.empty()) {
            int piv = v.back();
            if (basis[piv].empty()) {
                basis[piv] = std::move(v);
                ++rank;
                break;
            }
            xor_inplace(v, basis[piv], tmp);
        }
    }
    return rank;
}

static u64 C_mod(int a, int b, const std::vector<int> &spf, const std::vector<int> &prime_id, int num_primes,
                 u64 mod) {
    const int n = b - a + 1;
    const int r = rank_squarefree_vectors(a, b, spf, prime_id, num_primes);
    const int nullity = n - r;
    u64 total = mod_pow(2, (u64)nullity, mod);
    return (total + mod - 1) % mod;
}

int main() {
    static constexpr u64 MOD = 1000000007ULL;
    const int A = 1000000;
    const int B = 1234567;

    std::vector<int> spf = sieve_spf(B);
    std::vector<int> primes = list_primes_from_spf(spf);
    const int P = (int)primes.size();
    std::vector<int> prime_id(B + 1, -1);
    for (int i = 0; i < P; ++i) prime_id[primes[i]] = i;

#ifndef NDEBUG
    assert(C_mod(5, 10, spf, prime_id, P, MOD) == 3);
    assert(C_mod(40, 55, spf, prime_id, P, MOD) == 15);
    assert(C_mod(1000, 1234, spf, prime_id, P, MOD) == 975523611ULL);
#endif

    std::cout << C_mod(A, B, spf, prime_id, P, MOD) << "\n";
    return 0;
}

