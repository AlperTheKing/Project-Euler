#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using boost::multiprecision::cpp_int;

u64 mod_pow_u64(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return result;
}

int v_p_factorial(int n, int p) {
    int e = 0;
    while (n > 0) {
        n /= p;
        e += n;
    }
    return e;
}

int v2_u64(u64 x) {
    int v = 0;
    while ((x & 1ULL) == 0ULL) {
        x >>= 1ULL;
        ++v;
    }
    return v;
}

std::pair<std::vector<int>, std::vector<int>> build_spf_and_primes(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(p) * i;
            if (v > n || p > spf[i]) {
                break;
            }
            spf[static_cast<int>(v)] = p;
        }
    }

    return {std::move(spf), std::move(primes)};
}

std::vector<std::pair<int, int>> factorize_int(int x, const std::vector<int>& spf) {
    std::vector<std::pair<int, int>> fac;
    while (x > 1) {
        int p = spf[x];
        int e = 0;
        do {
            x /= p;
            ++e;
        } while (x > 1 && spf[x] == p);
        fac.push_back({p, e});
    }
    return fac;
}

u64 compute_R_mod(u64 p, int n, u64 out_mod) {
    const auto [spf, primes] = build_spf_and_primes(n);
    std::vector<int> max_exp(n + 1, 0);

    for (int q : primes) {
        if (q > n) {
            break;
        }

        const int e = v_p_factorial(n, q);

        if (q == 2) {
            int exp2 = 0;
            if (e >= 2) {
                if ((p & 3ULL) == 1ULL) {
                    const int u = v2_u64(p - 1ULL);
                    exp2 = (e <= u) ? 0 : (e - u);
                } else {
                    const int v = v2_u64(p + 1ULL);
                    exp2 = (e <= v) ? 1 : (e - v);
                }
            }
            if (exp2 > max_exp[2]) {
                max_exp[2] = exp2;
            }
            continue;
        }

        u64 t = static_cast<u64>(q - 1);
        const auto fac_qm1 = factorize_int(q - 1, spf);
        const u64 base_q = p % static_cast<u64>(q);

        for (const auto& [r, cnt] : fac_qm1) {
            for (int i = 0; i < cnt; ++i) {
                if (t % static_cast<u64>(r) == 0ULL &&
                    mod_pow_u64(base_q, t / static_cast<u64>(r), static_cast<u64>(q)) == 1ULL) {
                    t /= static_cast<u64>(r);
                } else {
                    break;
                }
            }
        }

        u64 tmp = t;
        for (const auto& [r, _] : fac_qm1) {
            int cnt = 0;
            while (tmp % static_cast<u64>(r) == 0ULL) {
                tmp /= static_cast<u64>(r);
                ++cnt;
            }
            if (cnt > max_exp[r]) {
                max_exp[r] = cnt;
            }
        }
        assert(tmp == 1ULL);

        int s = 1;
        if (e > 1) {
            const u64 mod_q2 = static_cast<u64>(q) * static_cast<u64>(q);
            if (mod_pow_u64(p % mod_q2, t, mod_q2) == 1ULL) {
                s = 2;
                if (e > 2) {
                    cpp_int base = p;
                    cpp_int exp = t;
                    cpp_int mod = cpp_int(q) * q * q;
                    for (int k = 3; k <= e; ++k) {
                        if (boost::multiprecision::powm(base, exp, mod) == 1) {
                            s = k;
                            mod *= q;
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        const int extra = e - s;
        if (extra > max_exp[q]) {
            max_exp[q] = extra;
        }
    }

    u64 ans = 1 % out_mod;
    for (int r = 2; r <= n; ++r) {
        if (max_exp[r] > 0) {
            ans = static_cast<u64>((static_cast<u128>(ans) * mod_pow_u64(static_cast<u64>(r), static_cast<u64>(max_exp[r]), out_mod)) % out_mod);
        }
    }

    return ans;
}

void run_validations() {
    assert(compute_R_mod(7ULL, 4, 1'000'000'007ULL) == 2ULL);
    assert(compute_R_mod(1'000'000'007ULL, 12, 1'000'000'007ULL) == 17'280ULL);
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kP = 1'000'000'007ULL;
    constexpr int kN = 10'000'000;
    std::cout << compute_R_mod(kP, kN, kP) << '\n';
    return 0;
}
