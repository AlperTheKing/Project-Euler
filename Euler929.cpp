#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr i64 kMod = 1'111'124'111LL;

struct NttMod {
    i64 mod;
    i64 primitive_root;
};

constexpr std::array<NttMod, 3> kNttMods = {
    NttMod{998'244'353LL, 3LL},
    NttMod{1'004'535'809LL, 3LL},
    NttMod{469'762'049LL, 3LL},
};

i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1LL) {
            r = static_cast<i64>((static_cast<u128>(r) * a) % mod);
        }
        a = static_cast<i64>((static_cast<u128>(a) * a) % mod);
        e >>= 1LL;
    }
    return r;
}

i64 mod_inv_prime(i64 a, i64 mod) {
    return mod_pow(a, mod - 2, mod);
}

void ntt(std::vector<i64>& a, i64 mod, i64 primitive_root, bool invert) {
    const int n = static_cast<int>(a.size());

    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        const i64 wlen = mod_pow(primitive_root, (mod - 1) / len, mod);
        const i64 wlen_inv = mod_inv_prime(wlen, mod);
        const i64 step = invert ? wlen_inv : wlen;

        for (int i = 0; i < n; i += len) {
            i64 w = 1;
            for (int j = 0; j < len / 2; ++j) {
                const i64 u = a[i + j];
                const i64 v = static_cast<i64>((static_cast<u128>(a[i + j + len / 2]) * w) % mod);
                i64 x = u + v;
                if (x >= mod) {
                    x -= mod;
                }
                i64 y = u - v;
                if (y < 0) {
                    y += mod;
                }
                a[i + j] = x;
                a[i + j + len / 2] = y;
                w = static_cast<i64>((static_cast<u128>(w) * step) % mod);
            }
        }
    }

    if (invert) {
        const i64 n_inv = mod_inv_prime(n, mod);
        for (i64& x : a) {
            x = static_cast<i64>((static_cast<u128>(x) * n_inv) % mod);
        }
    }
}

std::vector<i64> convolution_single_mod(const std::vector<i64>& a,
                                        const std::vector<i64>& b,
                                        i64 mod,
                                        i64 primitive_root) {
    if (a.empty() || b.empty()) {
        return {};
    }

    if (std::min(a.size(), b.size()) <= 64) {
        std::vector<i64> out(a.size() + b.size() - 1, 0);
        for (std::size_t i = 0; i < a.size(); ++i) {
            for (std::size_t j = 0; j < b.size(); ++j) {
                out[i + j] = (out[i + j] + static_cast<i64>((static_cast<u128>(a[i]) * b[j]) % mod)) % mod;
            }
        }
        return out;
    }

    int n = 1;
    while (n < static_cast<int>(a.size() + b.size() - 1)) {
        n <<= 1;
    }

    std::vector<i64> fa(n, 0), fb(n, 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        fa[i] = a[i] % mod;
    }
    for (std::size_t i = 0; i < b.size(); ++i) {
        fb[i] = b[i] % mod;
    }

    ntt(fa, mod, primitive_root, false);
    ntt(fb, mod, primitive_root, false);

    for (int i = 0; i < n; ++i) {
        fa[i] = static_cast<i64>((static_cast<u128>(fa[i]) * fb[i]) % mod);
    }

    ntt(fa, mod, primitive_root, true);
    fa.resize(a.size() + b.size() - 1);
    return fa;
}

std::vector<i64> convolution_mod(const std::vector<i64>& a, const std::vector<i64>& b, i64 target_mod) {
    if (a.empty() || b.empty()) {
        return {};
    }

    std::array<std::vector<i64>, 3> convs;
    for (int t = 0; t < 3; ++t) {
        convs[t] = convolution_single_mod(a, b, kNttMods[t].mod, kNttMods[t].primitive_root);
    }

    const i64 p1 = kNttMods[0].mod;
    const i64 p2 = kNttMods[1].mod;
    const i64 p3 = kNttMods[2].mod;

    const i64 inv_p1_mod_p2 = mod_inv_prime(p1 % p2, p2);
    const i64 p1_mod_p3 = p1 % p3;
    const i64 p1p2_mod_p3 = static_cast<i64>((static_cast<u128>(p1_mod_p3) * (p2 % p3)) % p3);
    const i64 inv_p1p2_mod_p3 = mod_inv_prime(p1p2_mod_p3, p3);

    const i64 p1_mod_m = p1 % target_mod;
    const i64 p1p2_mod_m = static_cast<i64>((static_cast<u128>(p1_mod_m) * (p2 % target_mod)) % target_mod);

    std::vector<i64> out(convs[0].size(), 0);
    for (std::size_t i = 0; i < out.size(); ++i) {
        const i64 r1 = convs[0][i];
        const i64 r2 = convs[1][i];
        const i64 r3 = convs[2][i];

        i64 t1 = r2 - r1;
        if (t1 < 0) {
            t1 += p2;
        }
        t1 = static_cast<i64>((static_cast<u128>(t1) * inv_p1_mod_p2) % p2);

        const i64 x12_mod_p3 = (r1 + static_cast<i64>((static_cast<u128>(p1_mod_p3) * t1) % p3)) % p3;

        i64 t2 = r3 - x12_mod_p3;
        if (t2 < 0) {
            t2 += p3;
        }
        t2 = static_cast<i64>((static_cast<u128>(t2) * inv_p1p2_mod_p3) % p3);

        i64 x = r1 % target_mod;
        x += static_cast<i64>((static_cast<u128>(p1_mod_m) * (t1 % target_mod)) % target_mod);
        if (x >= target_mod) {
            x -= target_mod;
        }
        x += static_cast<i64>((static_cast<u128>(p1p2_mod_m) * (t2 % target_mod)) % target_mod);
        x %= target_mod;
        out[i] = x;
    }

    return out;
}

std::vector<i64> poly_inv_one_minus(const std::vector<i64>& b, int n) {
    std::vector<i64> inv(1, 1);

    int len = 1;
    while (len < n) {
        const int m = std::min(2 * len, n);

        std::vector<i64> b_cut(m, 0);
        for (int i = 0; i < m; ++i) {
            b_cut[i] = b[i] % kMod;
        }

        std::vector<i64> prod = convolution_mod(inv, b_cut, kMod);
        prod.resize(m, 0);

        for (int i = 0; i < m; ++i) {
            prod[i] = (kMod - prod[i]) % kMod;
        }
        prod[0] = (prod[0] + 2) % kMod;

        std::vector<i64> next = convolution_mod(inv, prod, kMod);
        next.resize(m, 0);
        inv.swap(next);

        len = m;
    }

    inv.resize(n, 0);
    return inv;
}

std::vector<i64> build_A(int n) {
    std::vector<i64> fib(n + 1, 0);
    if (n >= 1) {
        fib[1] = 1;
    }
    for (int i = 2; i <= n; ++i) {
        fib[i] = fib[i - 1] + fib[i - 2];
        if (fib[i] >= kMod) {
            fib[i] -= kMod;
        }
    }

    std::vector<i64> a(n + 1, 0);

    for (int d = 1; d <= n; ++d) {
        i64 b = fib[d];
        if ((d & 1) == 0) {
            b = (kMod - b) % kMod;
        }

        for (int m = d; m <= n; m += d) {
            a[m] += b;
            if (a[m] >= kMod) {
                a[m] -= kMod;
            }
        }
    }

    return a;
}

i64 solve_fast(int n) {
    std::vector<i64> a = build_A(n);

    std::vector<i64> one_minus_a(n + 1, 0);
    one_minus_a[0] = 1;
    for (int i = 1; i <= n; ++i) {
        one_minus_a[i] = (kMod - a[i]) % kMod;
    }

    std::vector<i64> inv = poly_inv_one_minus(one_minus_a, n + 1);
    return inv[n] % kMod;
}

i64 solve_slow(int n) {
    std::vector<i64> a = build_A(n);
    std::vector<i64> u(n + 1, 0);
    u[0] = 1;

    for (int i = 1; i <= n; ++i) {
        i64 cur = 0;
        for (int k = 1; k <= i; ++k) {
            cur += static_cast<i64>((static_cast<u128>(a[k]) * u[i - k]) % kMod);
            if (cur >= kMod) {
                cur -= kMod;
            }
        }
        u[i] = cur;
    }

    return u[n];
}

void run_validations() {
    assert(solve_slow(5) == 10);

    for (int n = 1; n <= 80; ++n) {
        assert(solve_fast(n) == solve_slow(n));
    }
}

}  // namespace

int main() {
    run_validations();
    constexpr int kN = 100'000;
    std::cout << solve_fast(kN) << '\n';
    return 0;
}
