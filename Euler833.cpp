#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

using int64 = long long;
using u64 = unsigned long long;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr int64 kMod = 136101521;

struct PairInfo {
    int i_idx;
    int j_idx;
    u64 k_max;
    vector<int64> poly; // coefficients mod kMod
    int64 sum_all = 0;
    int64 sum_nf = 0;
};

u128 pow10_u128(int exp) {
    u128 res = 1;
    for (int i = 0; i < exp; ++i) res *= 10;
    return res;
}

bool mul_leq(u128 a, u128 b, u128 limit, u128 &prod) {
    if (a == 0 || b == 0) {
        prod = 0;
        return true;
    }
    if (a > limit / b) return false;
    prod = a * b;
    return true;
}

// Evaluate c = T_k * U_{i_idx}(x) * U_{j_idx}(x) with x=2k+1.
// Returns true if c <= limit.
bool eval_pair_leq(u64 k, int i_idx, int j_idx, u128 limit) {
    u128 kk = static_cast<u128>(k);
    u128 x = 2 * kk + 1;
    int max_idx = max(i_idx, j_idx);
    vector<u128> U(max_idx + 1, 0);
    U[0] = 1;
    if (max_idx >= 1) U[1] = 2 * x;

    u128 two_x = 2 * x;
    for (int m = 1; m < max_idx; ++m) {
        if (U[m] > (limit + U[m - 1]) / two_x) {
            U[m + 1] = limit + 1;
        } else {
            u128 term = two_x * U[m];
            u128 next = term - U[m - 1];
            U[m + 1] = (next > limit) ? (limit + 1) : next;
        }
    }

    u128 t = kk * (kk + 1) / 2;
    u128 tmp = 0;
    if (!mul_leq(t, U[i_idx], limit, tmp)) return false;
    if (!mul_leq(tmp, U[j_idx], limit, tmp)) return false;
    return true;
}

u128 eval_pair_value(u64 k, int i_idx, int j_idx) {
    u128 kk = static_cast<u128>(k);
    u128 x = 2 * kk + 1;
    int max_idx = max(i_idx, j_idx);
    vector<u128> U(max_idx + 1, 0);
    U[0] = 1;
    if (max_idx >= 1) U[1] = 2 * x;
    for (int m = 1; m < max_idx; ++m) {
        U[m + 1] = 2 * x * U[m] - U[m - 1];
    }
    u128 t = kk * (kk + 1) / 2;
    return t * U[i_idx] * U[j_idx];
}

u64 max_k_for_pair(u128 n, int i_idx, int j_idx) {
    u64 lo = 0;
    u64 hi = 1;
    while (eval_pair_leq(hi, i_idx, j_idx, n)) {
        if (hi > (numeric_limits<u64>::max() / 2)) break;
        hi *= 2;
    }
    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;
        if (eval_pair_leq(mid, i_idx, j_idx, n)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

vector<u64> generate_non_fundamentals(u64 k_max) {
    vector<u64> nf;
    if (k_max < 4) return nf;
    u64 k1_max = static_cast<u64>(sqrtl(static_cast<long double>(k_max) / 4.0L)) + 2;
    for (u64 k1 = 1; k1 <= k1_max; ++k1) {
        u128 x1 = 2 * static_cast<u128>(k1) + 1;
        u128 x_prev = 1;
        u128 x_cur = x1;
        while (true) {
            u128 x_next = 2 * x1 * x_cur - x_prev;
            u64 k_next = static_cast<u64>((x_next - 1) / 2);
            if (k_next > k_max) break;
            nf.push_back(k_next);
            x_prev = x_cur;
            x_cur = x_next;
        }
    }
    sort(nf.begin(), nf.end());
    nf.erase(unique(nf.begin(), nf.end()), nf.end());
    return nf;
}

int64 mod_pow(int64 base, int64 exp) {
    int64 res = 1 % kMod;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1) res = (static_cast<i128>(res) * base) % kMod;
        base = (static_cast<i128>(base) * base) % kMod;
        exp >>= 1;
    }
    return res;
}

int64 mod_inv(int64 a) {
    return mod_pow(a, kMod - 2);
}

using Poly = vector<int64>;

Poly poly_add(const Poly &a, const Poly &b) {
    size_t n = max(a.size(), b.size());
    Poly c(n, 0);
    for (size_t i = 0; i < n; ++i) {
        int64 va = (i < a.size()) ? a[i] : 0;
        int64 vb = (i < b.size()) ? b[i] : 0;
        int64 v = va + vb;
        if (v >= kMod) v -= kMod;
        c[i] = v;
    }
    return c;
}

Poly poly_sub(const Poly &a, const Poly &b) {
    size_t n = max(a.size(), b.size());
    Poly c(n, 0);
    for (size_t i = 0; i < n; ++i) {
        int64 va = (i < a.size()) ? a[i] : 0;
        int64 vb = (i < b.size()) ? b[i] : 0;
        int64 v = va - vb;
        if (v < 0) v += kMod;
        c[i] = v;
    }
    return c;
}

Poly poly_mul(const Poly &a, const Poly &b) {
    Poly c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            c[i + j] = (c[i + j] + static_cast<i128>(a[i]) * b[j]) % kMod;
        }
    }
    return c;
}

Poly poly_scale(const Poly &a, int64 s) {
    Poly c(a.size(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
        c[i] = static_cast<i128>(a[i]) * s % kMod;
    }
    return c;
}

vector<Poly> build_U_polys(int max_idx) {
    Poly t = {1, 2}; // 2k+1
    vector<Poly> U(max_idx + 1);
    U[0] = {1};
    if (max_idx >= 1) {
        U[1] = poly_scale(t, 2);
    }
    for (int m = 1; m < max_idx; ++m) {
        Poly temp = poly_mul(t, U[m]);
        temp = poly_scale(temp, 2);
        U[m + 1] = poly_sub(temp, U[m - 1]);
    }
    return U;
}

struct SumPowers {
    vector<vector<int64>> y;
    vector<int64> fact;
    vector<int64> invfact;
};

SumPowers precompute_sum_pows(int max_p) {
    SumPowers sp;
    sp.y.resize(max_p + 1);
    for (int p = 0; p <= max_p; ++p) {
        int m = p + 1;
        sp.y[p].assign(m + 1, 0);
        int64 s = 0;
        for (int i = 1; i <= m; ++i) {
            s += mod_pow(i, p);
            s %= kMod;
            sp.y[p][i] = s;
        }
    }
    sp.fact.assign(max_p + 2, 1);
    sp.invfact.assign(max_p + 2, 1);
    for (int i = 1; i < static_cast<int>(sp.fact.size()); ++i) {
        sp.fact[i] = static_cast<i128>(sp.fact[i - 1]) * i % kMod;
    }
    sp.invfact.back() = mod_inv(sp.fact.back());
    for (int i = static_cast<int>(sp.fact.size()) - 1; i >= 1; --i) {
        sp.invfact[i - 1] = static_cast<i128>(sp.invfact[i]) * i % kMod;
    }
    return sp;
}

int64 sum_pows(int p, u64 n, const SumPowers &sp) {
    if (n <= static_cast<u64>(p + 1)) {
        return sp.y[p][static_cast<size_t>(n)];
    }
    int m = p + 1;
    vector<int64> pre(m + 2, 1), suf(m + 2, 1);
    int64 nmod = static_cast<int64>(n % kMod);
    for (int i = 0; i <= m; ++i) {
        int64 term = nmod - i;
        if (term < 0) term += kMod;
        pre[i + 1] = static_cast<i128>(pre[i]) * term % kMod;
    }
    for (int i = m; i >= 0; --i) {
        int64 term = nmod - i;
        if (term < 0) term += kMod;
        suf[i] = static_cast<i128>(suf[i + 1]) * term % kMod;
    }
    int64 res = 0;
    for (int i = 0; i <= m; ++i) {
        int64 num = static_cast<i128>(pre[i]) * suf[i + 1] % kMod;
        int64 denom_inv = static_cast<i128>(sp.invfact[i]) * sp.invfact[m - i] % kMod;
        int64 term = static_cast<i128>(sp.y[p][i]) * num % kMod;
        term = static_cast<i128>(term) * denom_inv % kMod;
        if ((m - i) & 1) {
            if (term != 0) term = kMod - term;
        }
        res += term;
        if (res >= kMod) res -= kMod;
    }
    return res;
}

int64 sum_poly(const Poly &coeffs, u64 n, const SumPowers &sp) {
    int64 res = 0;
    for (size_t p = 0; p < coeffs.size(); ++p) {
        if (coeffs[p] == 0) continue;
        int64 spow = sum_pows(static_cast<int>(p), n, sp);
        res = (res + static_cast<i128>(coeffs[p]) * spow) % kMod;
    }
    return res;
}

vector<PairInfo> build_pairs(u128 n) {
    vector<u128> U3;
    U3.push_back(1);
    U3.push_back(6);
    while (true) {
        u128 next = 2 * static_cast<u128>(3) * U3[U3.size() - 1] - U3[U3.size() - 2];
        U3.push_back(next);
        if (next > n) break;
    }
    int max_idx = static_cast<int>(U3.size()) - 1;

    vector<PairInfo> pairs;
    for (int i = 0; i <= max_idx; ++i) {
        for (int j = i + 1; j <= max_idx; ++j) {
            u128 prod = 0;
            if (!mul_leq(U3[i], U3[j], n, prod)) continue;
            PairInfo info;
            info.i_idx = i;
            info.j_idx = j;
            pairs.push_back(info);
        }
    }
    return pairs;
}

u64 compute_exact_small(u64 n) {
    u128 n128 = static_cast<u128>(n);
    vector<PairInfo> pairs = build_pairs(n128);
    u64 k_max = 0;
    for (auto &p : pairs) {
        p.k_max = max_k_for_pair(n128, p.i_idx, p.j_idx);
        k_max = max(k_max, p.k_max);
    }
    vector<u64> nf = generate_non_fundamentals(k_max);
    vector<char> is_nf(k_max + 1, 0);
    for (u64 v : nf) if (v <= k_max) is_nf[v] = 1;

    u128 total = 0;
    for (const auto &p : pairs) {
        for (u64 k = 1; k <= p.k_max; ++k) {
            if (is_nf[k]) continue;
            u128 c = eval_pair_value(k, p.i_idx, p.j_idx);
            if (c > n128) break;
            total += c;
        }
    }
    return static_cast<u64>(total);
}

int64 compute_mod(u128 n) {
    vector<PairInfo> pairs = build_pairs(n);
    int max_idx = 0;
    u64 k_max = 0;
    for (auto &p : pairs) {
        p.k_max = max_k_for_pair(n, p.i_idx, p.j_idx);
        k_max = max(k_max, p.k_max);
        max_idx = max(max_idx, max(p.i_idx, p.j_idx));
    }

    vector<u64> nf = generate_non_fundamentals(k_max);

    const int64 inv2 = mod_inv(2);
    vector<Poly> U_polys = build_U_polys(max_idx);
    Poly t_poly = {0, inv2, inv2};
    int deg_max = 0;
    for (auto &p : pairs) {
        Poly r = poly_mul(U_polys[p.i_idx], U_polys[p.j_idx]);
        p.poly = poly_mul(r, t_poly);
        deg_max = max(deg_max, static_cast<int>(p.poly.size()) - 1);
    }

    SumPowers sp = precompute_sum_pows(deg_max);
    for (auto &p : pairs) {
        p.sum_all = sum_poly(p.poly, p.k_max, sp);
    }

    // Parallel non-fundamental contributions.
    const unsigned int hw = thread::hardware_concurrency();
    const unsigned int threads = (hw == 0) ? 4u : min(hw, 8u);
    vector<vector<int64>> local_sums(threads, vector<int64>(pairs.size(), 0));
    vector<thread> workers;

    auto worker = [&](unsigned int tid, size_t start, size_t end) {
        vector<int64> U(max_idx + 1, 0);
        for (size_t idx = start; idx < end; ++idx) {
            u64 k = nf[idx];
            u64 k_mod = k % kMod;
            int64 x_mod = static_cast<int64>((2 * (k_mod % kMod) + 1) % kMod);
            U[0] = 1;
            if (max_idx >= 1) U[1] = (2LL * x_mod) % kMod;
            for (int m = 1; m < max_idx; ++m) {
                int64 next = (2LL * x_mod % kMod * U[m] - U[m - 1]) % kMod;
                if (next < 0) next += kMod;
                U[m + 1] = next;
            }
            int64 t = static_cast<int64>(k_mod) * static_cast<int64>((k_mod + 1) % kMod) % kMod;
            t = static_cast<int64>(static_cast<i128>(t) * inv2 % kMod);

            auto &sum_nf = local_sums[tid];
            for (size_t p = 0; p < pairs.size(); ++p) {
                if (k > pairs[p].k_max) continue;
                int64 val = static_cast<int64>(static_cast<i128>(t) * U[pairs[p].i_idx] % kMod);
                val = static_cast<int64>(static_cast<i128>(val) * U[pairs[p].j_idx] % kMod);
                sum_nf[p] += val;
                if (sum_nf[p] >= kMod) sum_nf[p] -= kMod;
            }
        }
    };

    size_t chunk = (nf.size() + threads - 1) / threads;
    for (unsigned int t = 0; t < threads; ++t) {
        size_t start = t * chunk;
        size_t end = min(nf.size(), start + chunk);
        workers.emplace_back(worker, t, start, end);
    }
    for (auto &th : workers) th.join();

    for (size_t p = 0; p < pairs.size(); ++p) {
        int64 total = 0;
        for (unsigned int t = 0; t < threads; ++t) {
            total += local_sums[t][p];
            if (total >= kMod) total -= kMod;
        }
        pairs[p].sum_nf = total;
    }

    int64 ans = 0;
    for (const auto &p : pairs) {
        int64 term = p.sum_all - p.sum_nf;
        if (term < 0) term += kMod;
        ans += term;
        if (ans >= kMod) ans -= kMod;
    }
    return ans;
}

} // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Validation checks.
    if (compute_exact_small(100000) != 1479802ULL) {
        cerr << "Validation failed: S(1e5) mismatch\n";
        return 1;
    }
    if (compute_exact_small(1000000000ULL) != 241614948794ULL) {
        cerr << "Validation failed: S(1e9) mismatch\n";
        return 1;
    }

    u128 n = pow10_u128(35);
    int64 ans = compute_mod(n);
    cout << ans << "\n";
    return 0;
}
