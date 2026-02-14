#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr u64 L = 10000000000000000ULL;
constexpr int K = 10000000;
constexpr u64 MOD = 16ULL * 9ULL * 5ULL * 7ULL * 11ULL * 13ULL;
constexpr u64 START = (L / MOD) * MOD;
constexpr std::array<u64, 9> D = {1, 6, 1, 4, 315, 2, 1, 24, 1};

u64 egcd_inv(u64 a, u64 mod) {
    i64 t = 0, nt = 1;
    i64 r = static_cast<i64>(mod), nr = static_cast<i64>(a % mod);
    while (nr != 0) {
        i64 q = r / nr;
        i64 tt = t - q * nt;
        t = nt;
        nt = tt;
        i64 rr = r - q * nr;
        r = nr;
        nr = rr;
    }
    if (r != 1) return 0;
    if (t < 0) t += static_cast<i64>(mod);
    return static_cast<u64>(t);
}

u64 crt_two(u64 m1, u64 r1, u64 m2, u64 r2) {
    const u64 inv = egcd_inv(m1, m2);
    const u64 delta = (r2 >= r1 % m2) ? (r2 - (r1 % m2)) : (r2 + m2 - (r1 % m2));
    const u64 t = (static_cast<__uint128_t>(delta) * inv) % m2;
    return r1 + m1 * t;
}

std::vector<u64> build_residues() {
    std::vector<u64> rems = {1, 7};
    u64 prod = 16;

    for (u64& r : rems) r = crt_two(prod, r, 9, 5);
    prod *= 9;

    for (u64& r : rems) r = crt_two(prod, r, 5, 1);
    prod *= 5;

    for (u64& r : rems) r = crt_two(prod, r, 7, 3);
    prod *= 7;

    {
        std::vector<u64> nxt;
        nxt.reserve(rems.size() * 2);
        for (u64 r : rems) {
            nxt.push_back(crt_two(prod, r, 11, 1));
            nxt.push_back(crt_two(prod, r, 11, 2));
        }
        rems.swap(nxt);
    }
    prod *= 11;

    {
        std::vector<u64> nxt;
        nxt.reserve(rems.size() * 4);
        for (u64 r : rems) {
            nxt.push_back(crt_two(prod, r, 13, 1));
            nxt.push_back(crt_two(prod, r, 13, 2));
            nxt.push_back(crt_two(prod, r, 13, 3));
            nxt.push_back(crt_two(prod, r, 13, 4));
        }
        rems.swap(nxt);
    }

    std::sort(rems.begin(), rems.end());
    rems.erase(std::unique(rems.begin(), rems.end()), rems.end());
    return rems;
}

struct PrimeInv {
    u32 p;
    u32 inv;
};

std::vector<PrimeInv> build_primes_and_inverses() {
    constexpr int LIM = 100000000;
    std::vector<u8> is_prime(LIM, 1);
    is_prime[0] = 0;
    is_prime[1] = 0;
    for (int i = 2; static_cast<i64>(i) * i < LIM; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j < LIM; j += i) is_prime[j] = 0;
    }

    std::vector<PrimeInv> out;
    out.reserve(6000000);
    for (int p = 17; p < LIM; ++p) {
        if (!is_prime[p]) continue;
        if (p == 2 || p == 3 || p == 5 || p == 7 || p == 11 || p == 13) continue;
        const u64 inv = egcd_inv(MOD % static_cast<u64>(p), static_cast<u64>(p));
        out.push_back({static_cast<u32>(p), static_cast<u32>(inv)});
    }
    return out;
}

u64 solve() {
    const auto rems = build_residues();
    const auto pinv = build_primes_and_inverses();

    std::vector<u8> flag(K, 255);
    u64 best_sum = 0;

    for (u8 rid = 0; rid < rems.size(); ++rid) {
        const u64 rem = rems[rid];

        for (const auto [p, inv] : pinv) {
            u64 base = ((START + rem) % p);
            base = (static_cast<__uint128_t>(base) * inv) % p;
            for (int i = 0; i < 9; ++i) {
                for (u64 k = base; k < K; k += p) flag[static_cast<size_t>(k)] = rid;
                base += inv;
                if (base >= p) base -= p;
            }
        }

        for (u64 i = 0; i < K; ++i) {
            if (flag[static_cast<size_t>(i)] == rid) continue;
            const u64 cand = START + rem - MOD * i;
            u64 sum = 0;
            for (int j = 0; j < 9; ++j) sum += (cand + static_cast<u64>(j)) / D[j];
            if (sum > best_sum) best_sum = sum;
        }
    }

    return best_sum;
}

u64 lpf_small(u64 n) {
    u64 best = 1;
    for (u64 p = 2; p * p <= n; ++p) {
        while (n % p == 0) {
            best = p;
            n /= p;
        }
    }
    if (n > 1) best = n;
    return best;
}

bool checks() {
    if (lpf_small(100) != 5ULL) return false;
    u64 g100 = 0;
    for (u64 n = 100; n <= 108; ++n) g100 += lpf_small(n);
    if (g100 != 409ULL) return false;
    u64 h100 = 0;
    for (u64 x = 2; x <= 100; ++x) {
        u64 s = 0;
        for (u64 n = x; n <= x + 8; ++n) s += lpf_small(n);
        h100 = std::max(h100, s);
    }
    if (h100 != 417ULL) return false;
    return true;
}

int main() {
    if (!checks()) return 1;
    std::cout << solve() << '\n';
}
