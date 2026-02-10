#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

// Project Euler 598: Split Divisibilities
//
// For n = 100!, write n = \prod p_i^{e_i}. A divisor a corresponds to exponents x_i (0..e_i),
// and b = n/a has exponents e_i - x_i.
// The divisor counts are:
//   tau(a) = \prod (x_i + 1),   tau(b) = \prod (e_i - x_i + 1).
// We count exponent vectors with tau(a) = tau(b); each non-fixed solution pairs with its
// complement, so C(n) = (solutions + fixed)/2 where fixed exists iff all e_i are even.
//
// For 100!, primes > 7 have exponent <= 9, so (x+1) and (e-x+1) are <= 10 and factor only
// into primes 2,3,5,7. Therefore all primes > 7 must cancel internally among the big primes
// 2,3,5,7. We enumerate the 2*3*5*7 part by brute force (about 2 million choices), filter
// those with zero contribution on primes > 7, and record the resulting exponent differences
// for primes 2,3,5,7. The remaining primes are handled by a small DP on those 4 exponents.

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    std::vector<int> ps;
    for (int i = 2; i <= n; ++i)
        if (is_prime[i]) ps.push_back(i);
    return ps;
}

static int exp_in_fact(int n, int p) {
    int e = 0;
    while (n) {
        n /= p;
        e += n;
    }
    return e;
}

static u64 pack4(int d2, int d3, int d5, int d7) {
    // Pack signed values into 4x16-bit lanes.
    static constexpr int OFF = 128;
    const u64 a = (u64)(std::uint16_t)(d2 + OFF);
    const u64 b = (u64)(std::uint16_t)(d3 + OFF);
    const u64 c = (u64)(std::uint16_t)(d5 + OFF);
    const u64 d = (u64)(std::uint16_t)(d7 + OFF);
    return a | (b << 16) | (c << 32) | (d << 48);
}

static void unpack4(u64 key, int& d2, int& d3, int& d5, int& d7) {
    static constexpr int OFF = 128;
    d2 = (int)(std::uint16_t)(key & 0xFFFFULL) - OFF;
    d3 = (int)(std::uint16_t)((key >> 16) & 0xFFFFULL) - OFF;
    d5 = (int)(std::uint16_t)((key >> 32) & 0xFFFFULL) - OFF;
    d7 = (int)(std::uint16_t)((key >> 48) & 0xFFFFULL) - OFF;
}

static u128 C_factorial(int n) {
    const auto ps = primes_up_to(n);

    // Exponents in n!.
    std::vector<int> exps;
    exps.reserve(ps.size());
    for (int p : ps) exps.push_back(exp_in_fact(n, p));

    // Fixed point exists iff all exponents are even.
    bool all_even = true;
    for (int e : exps)
        if (e & 1) {
            all_even = false;
            break;
        }

    // Big primes handled by brute force: p in {2,3,5,7} (assumed present for n>=7).
    const std::array<int, 4> bigP{2, 3, 5, 7};

    // For n<=100, all primes > 7 have exponent <= 9.
    for (size_t i = 0; i < ps.size(); ++i) {
        if (ps[i] > 7) assert(exps[i] <= 9);
    }

    int max_val = 1;
    for (size_t i = 0; i < ps.size(); ++i) {
        max_val = std::max(max_val, exps[i] + 1);
    }

    // Prime basis up to max_val.
    const auto basis = primes_up_to(max_val);
    const int P = (int)basis.size();
    std::unordered_map<int, int> p2idx;
    p2idx.reserve(basis.size());
    for (int i = 0; i < P; ++i) p2idx[basis[i]] = i;

    // Factorization exponent vectors for 1..max_val.
    std::vector<std::vector<std::int8_t>> fac((size_t)max_val + 1, std::vector<std::int8_t>((size_t)P, 0));
    for (int x = 2; x <= max_val; ++x) {
        int m = x;
        for (int i = 0; i < P; ++i) {
            const int p = basis[i];
            if (p * p > m) break;
            while (m % p == 0) {
                ++fac[(size_t)x][(size_t)i];
                m /= p;
            }
        }
        if (m > 1) {
            const int i = p2idx[m];
            ++fac[(size_t)x][(size_t)i];
        }
    }

    const int idx2 = p2idx[2];
    const int idx3 = p2idx[3];
    const int idx5 = p2idx[5];
    const int idx7 = p2idx[7];

    std::vector<int> big_indices;
    for (int i = 0; i < P; ++i) {
        if (basis[i] > 7) big_indices.push_back(i);
    }

    // Build diffs for each big prime variable.
    std::vector<std::vector<std::vector<std::int8_t>>> diffs;
    diffs.reserve(4);
    for (int bp : bigP) {
        if (bp > n) {
            diffs.push_back({});
            continue;
        }
        int e = exp_in_fact(n, bp);
        const int E = e + 2;
        std::vector<std::vector<std::int8_t>> v;
        v.reserve((size_t)(E - 1));
        for (int y = 1; y <= e + 1; ++y) {
            const int z = E - y;
            std::vector<std::int8_t> dv((size_t)P, 0);
            for (int i = 0; i < P; ++i) dv[(size_t)i] = (std::int8_t)(fac[(size_t)y][(size_t)i] - fac[(size_t)z][(size_t)i]);
            v.push_back(std::move(dv));
        }
        diffs.push_back(std::move(v));
    }

    // Pair sums (0,1) and (2,3).
    std::vector<std::vector<std::int8_t>> pair01;
    std::vector<std::vector<std::int8_t>> pair23;

    pair01.reserve(diffs[0].size() * diffs[1].size());
    for (const auto& a : diffs[0]) {
        for (const auto& b : diffs[1]) {
            std::vector<std::int8_t> s((size_t)P, 0);
            for (int i = 0; i < P; ++i) s[(size_t)i] = (std::int8_t)(a[(size_t)i] + b[(size_t)i]);
            pair01.push_back(std::move(s));
        }
    }

    pair23.reserve(diffs[2].size() * diffs[3].size());
    for (const auto& a : diffs[2]) {
        for (const auto& b : diffs[3]) {
            std::vector<std::int8_t> s((size_t)P, 0);
            for (int i = 0; i < P; ++i) s[(size_t)i] = (std::int8_t)(a[(size_t)i] + b[(size_t)i]);
            pair23.push_back(std::move(s));
        }
    }

    // map1: big group assignments with all primes>7 canceled, keyed by (d2,d3,d5,d7).
    std::unordered_map<u64, u64> map1;
    map1.reserve(10000);

    for (const auto& a : pair01) {
        for (const auto& b : pair23) {
            bool ok = true;
            for (int i : big_indices) {
                if ((int)a[(size_t)i] + (int)b[(size_t)i] != 0) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            const int d2 = (int)a[(size_t)idx2] + (int)b[(size_t)idx2];
            const int d3 = (int)a[(size_t)idx3] + (int)b[(size_t)idx3];
            const int d5 = (int)a[(size_t)idx5] + (int)b[(size_t)idx5];
            const int d7 = (int)a[(size_t)idx7] + (int)b[(size_t)idx7];
            ++map1[pack4(d2, d3, d5, d7)];
        }
    }

    // map2: DP over primes > 7, but only tracking exponents for 2,3,5,7.
    std::unordered_map<u64, u64> map2;
    map2.reserve(30000);
    map2[pack4(0, 0, 0, 0)] = 1;

    for (size_t i = 0; i < ps.size(); ++i) {
        const int p = ps[i];
        if (p <= 7) continue;
        const int e = exps[i];
        const int E = e + 2;

        std::vector<std::array<int, 4>> deltas;
        deltas.reserve((size_t)(E - 1));
        for (int y = 1; y <= e + 1; ++y) {
            const int z = E - y;
            deltas.push_back({(int)fac[(size_t)y][(size_t)idx2] - (int)fac[(size_t)z][(size_t)idx2],
                              (int)fac[(size_t)y][(size_t)idx3] - (int)fac[(size_t)z][(size_t)idx3],
                              (int)fac[(size_t)y][(size_t)idx5] - (int)fac[(size_t)z][(size_t)idx5],
                              (int)fac[(size_t)y][(size_t)idx7] - (int)fac[(size_t)z][(size_t)idx7]});
        }

        std::unordered_map<u64, u64> nxt;
        nxt.reserve(map2.size() * deltas.size());
        for (const auto& kv : map2) {
            int d2, d3, d5, d7;
            unpack4(kv.first, d2, d3, d5, d7);
            const u64 cnt = kv.second;
            for (const auto& dd : deltas) {
                const u64 nk = pack4(d2 + dd[0], d3 + dd[1], d5 + dd[2], d7 + dd[3]);
                nxt[nk] += cnt;
            }
        }
        map2.swap(nxt);
    }

    u128 M = 0;
    for (const auto& kv : map1) {
        int d2, d3, d5, d7;
        unpack4(kv.first, d2, d3, d5, d7);
        const u64 need = pack4(-d2, -d3, -d5, -d7);
        auto it = map2.find(need);
        if (it == map2.end()) continue;
        M += (u128)kv.second * (u128)it->second;
    }

    const u128 fixed = all_even ? 1 : 0;
    return (M + fixed) / 2;
}

static void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }
    char buf[64];
    int n = 0;
    while (x > 0) {
        buf[n++] = (char)('0' + (u64)(x % 10));
        x /= 10;
    }
    while (n--) std::cout << buf[n];
}

int main() {
    // Validation from the statement.
    if (C_factorial(10) != 3) {
        std::cerr << "Validation failed: C(10!)\n";
        return 1;
    }

    const u128 ans = C_factorial(100);
    print_u128(ans);
    std::cout << "\n";
    return 0;
}
