#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) r = static_cast<u64>((static_cast<u128>(r) * a) % mod);
        a = static_cast<u64>((static_cast<u128>(a) * a) % mod);
        e >>= 1ULL;
    }
    return r;
}

static u64 tonelli_shanks(u64 n, u64 p) {
    if (n == 0) return 0;
    if (p == 2) return n;
    if (mod_pow(n, (p - 1) / 2, p) != 1) return 0;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);

    u64 q = p - 1;
    u64 s = 0;
    while ((q & 1ULL) == 0) {
        q >>= 1ULL;
        ++s;
    }

    u64 z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) ++z;

    u64 c = mod_pow(z, q, p);
    u64 x = mod_pow(n, (q + 1) / 2, p);
    u64 t = mod_pow(n, q, p);
    u64 m = s;

    while (t != 1) {
        u64 i = 1;
        u64 tt = static_cast<u64>((static_cast<u128>(t) * t) % p);
        while (tt != 1) {
            tt = static_cast<u64>((static_cast<u128>(tt) * tt) % p);
            ++i;
        }
        u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        x = static_cast<u64>((static_cast<u128>(x) * b) % p);
        c = static_cast<u64>((static_cast<u128>(b) * b) % p);
        t = static_cast<u64>((static_cast<u128>(t) * c) % p);
        m = i;
    }
    return x;
}

static u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) > 0 && (r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        s.push_back(static_cast<char>('0' + (x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u128 solve_F(u64 N) {
    std::vector<bool> is_prime(N + 1, true);
    is_prime[0] = false;
    if (N >= 1) is_prime[1] = false;
    for (u64 i = 2; i * i <= N; ++i) {
        if (!is_prime[i]) continue;
        for (u64 j = i * i; j <= N; j += i) is_prime[j] = false;
    }

    std::vector<int> root_odd(N + 1, -1);
    std::vector<int> vals;
    vals.reserve(70000);
    vals.push_back(1);
    if (N >= 2) vals.push_back(2);

    for (u64 p = 3; p <= N; p += 2) {
        if (!is_prime[p] || p % 4 != 1) continue;
        u64 r = tonelli_shanks(p - 1, p);
        u64 m = p;
        while (m <= N) {
            if (root_odd[m] == -1) root_odd[m] = static_cast<int>(r);
            vals.push_back(static_cast<int>(m));
            if (2 * m <= N) vals.push_back(static_cast<int>(2 * m));

            if (m > N / p) break;
            u64 f_div_m = static_cast<u64>(((static_cast<u128>(r) * r + 1) / m) % p);
            u64 inv = mod_pow((2 * (r % p)) % p, p - 2, p);
            u64 t = (p - static_cast<u64>((static_cast<u128>(f_div_m) * inv) % p)) % p;
            r += t * m;
            m *= p;
            r %= m;
        }
    }

    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());

    std::vector<int> id(N + 1, -1);
    for (int i = 0; i < static_cast<int>(vals.size()); ++i) id[vals[i]] = i;

    std::vector<char> allowed(N + 1, 0);
    for (int v : vals) allowed[v] = 1;

    std::vector<std::vector<int>> adj(vals.size());
    for (int va : vals) {
        std::vector<u64> roots;
        if (va == 1) {
            roots.push_back(0);
        } else if (va == 2) {
            roots.push_back(1);
        } else if ((va & 1) == 0) {
            u64 m = static_cast<u64>(va / 2);
            int rr = root_odd[m];
            if (rr == -1) continue;
            u64 r1 = static_cast<u64>(rr);
            u64 r2 = m - r1;
            u64 mod = 2 * m;
            if ((r1 & 1ULL) == 0) r1 += m;
            if ((r2 & 1ULL) == 0) r2 += m;
            r1 %= mod;
            r2 %= mod;
            roots.push_back(r1);
            if (r2 != r1) roots.push_back(r2);
        } else {
            int rr = root_odd[va];
            if (rr == -1) continue;
            u64 r1 = static_cast<u64>(rr);
            u64 r2 = static_cast<u64>(va) - r1;
            roots.push_back(r1);
            if (r2 != r1) roots.push_back(r2);
        }

        u64 limit = isqrt_u64(static_cast<u64>(static_cast<u128>(va) * N - 1));
        int ia = id[va];
        for (u64 r : roots) {
            u64 x = r;
            if (x == 0) x += va;
            for (; x <= limit; x += va) {
                u64 b = (x * x + 1) / static_cast<u64>(va);
                if (b > N || b <= static_cast<u64>(va)) continue;
                if (!allowed[b]) continue;
                int ib = id[b];
                if (ib < 0 || ib == ia) continue;
                adj[ia].push_back(ib);
                adj[ib].push_back(ia);
            }
        }
    }

    for (auto& v : adj) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    std::vector<char> vis(vals.size(), 0);
    u128 total_oriented = 0;

    std::function<void(int, int, int, u128)> dfs = [&](int start, int u, int depth, u128 sumv) {
        for (int v : adj[u]) {
            if (v < start) continue;
            if (v == start) {
                if (depth >= 3) total_oriented += sumv;
            } else if (!vis[v]) {
                vis[v] = 1;
                dfs(start, v, depth + 1, sumv + static_cast<u128>(vals[v]));
                vis[v] = 0;
            }
        }
    };

    for (int s = 0; s < static_cast<int>(vals.size()); ++s) {
        if (adj[s].size() < 2) continue;
        vis[s] = 1;
        dfs(s, s, 1, static_cast<u128>(vals[s]));
        vis[s] = 0;
    }

    return total_oriented / 2;
}

int main() {
    assert(solve_F(20) == 258);
    assert(solve_F(100) == 538768);
    std::cout << to_string_u128(solve_F(1'000'000)) << '\n';
    return 0;
}
