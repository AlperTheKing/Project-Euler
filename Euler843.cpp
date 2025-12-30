#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

namespace {

struct Poly {
    uint64_t lo = 0;
    uint64_t hi = 0;

    bool is_zero() const {
        return lo == 0 && hi == 0;
    }
    bool is_one() const {
        return lo == 1 && hi == 0;
    }
};

Poly operator^(const Poly& a, const Poly& b) {
    return {a.lo ^ b.lo, a.hi ^ b.hi};
}

Poly& operator^=(Poly& a, const Poly& b) {
    a.lo ^= b.lo;
    a.hi ^= b.hi;
    return a;
}

bool poly_equal(const Poly& a, const Poly& b) {
    return a.lo == b.lo && a.hi == b.hi;
}

int poly_degree(const Poly& p) {
    if (p.hi != 0) {
        return 64 + 63 - __builtin_clzll(p.hi);
    }
    if (p.lo != 0) {
        return 63 - __builtin_clzll(p.lo);
    }
    return -1;
}

bool poly_get_bit(const Poly& p, int pos) {
    if (pos < 64) return (p.lo >> pos) & 1ULL;
    return (p.hi >> (pos - 64)) & 1ULL;
}

void poly_toggle_bit(Poly& p, int pos) {
    if (pos < 64) {
        p.lo ^= 1ULL << pos;
    } else {
        p.hi ^= 1ULL << (pos - 64);
    }
}

Poly poly_shift_left(const Poly& p, int k) {
    unsigned __int128 v = (static_cast<unsigned __int128>(p.hi) << 64) | p.lo;
    v <<= k;
    return {static_cast<uint64_t>(v), static_cast<uint64_t>(v >> 64)};
}

Poly poly_shift_right1(const Poly& p) {
    Poly r;
    r.lo = (p.lo >> 1) | (p.hi << 63);
    r.hi = p.hi >> 1;
    return r;
}

Poly poly_shift_left1(const Poly& p) {
    Poly r;
    r.hi = (p.hi << 1) | (p.lo >> 63);
    r.lo = p.lo << 1;
    return r;
}

Poly poly_mod(Poly a, const Poly& mod) {
    int deg_mod = poly_degree(mod);
    if (deg_mod < 0) return a;
    for (;;) {
        int deg_a = poly_degree(a);
        if (deg_a < deg_mod) break;
        int shift = deg_a - deg_mod;
        a ^= poly_shift_left(mod, shift);
    }
    return a;
}

Poly poly_div(Poly a, const Poly& mod) {
    Poly q;
    int deg_mod = poly_degree(mod);
    for (;;) {
        int deg_a = poly_degree(a);
        if (deg_a < deg_mod) break;
        int shift = deg_a - deg_mod;
        poly_toggle_bit(q, shift);
        a ^= poly_shift_left(mod, shift);
    }
    return q;
}

Poly poly_gcd(Poly a, Poly b) {
    while (!b.is_zero()) {
        Poly r = poly_mod(a, b);
        a = b;
        b = r;
    }
    return a;
}

Poly poly_mul_mod(Poly a, Poly b, const Poly& mod) {
    if (a.is_zero() || b.is_zero()) return {0, 0};
    a = poly_mod(a, mod);
    b = poly_mod(b, mod);
    int deg_mod = poly_degree(mod);
    Poly res;
    while (!b.is_zero()) {
        if (b.lo & 1ULL) res ^= a;
        b = poly_shift_right1(b);
        a = poly_shift_left1(a);
        if (poly_get_bit(a, deg_mod)) a ^= mod;
    }
    return res;
}

Poly poly_pow_mod(Poly base, uint64_t exp, const Poly& mod) {
    Poly res{1, 0};
    while (exp > 0) {
        if (exp & 1ULL) res = poly_mul_mod(res, base, mod);
        base = poly_mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return res;
}

vector<uint32_t> sieve_primes(int limit) {
    vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= limit; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= limit; j += i) is_prime[j] = false;
    }
    vector<uint32_t> primes;
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[i]) primes.push_back(static_cast<uint32_t>(i));
    }
    return primes;
}

vector<uint64_t> factorize_uint64(uint64_t n, const vector<uint32_t>& primes) {
    vector<uint64_t> factors;
    uint64_t temp = n;
    for (uint32_t p : primes) {
        uint64_t pp = static_cast<uint64_t>(p);
        if (pp * pp > temp) break;
        if (temp % pp == 0) {
            factors.push_back(pp);
            while (temp % pp == 0) temp /= pp;
        }
    }
    if (temp > 1) factors.push_back(temp);
    return factors;
}

vector<Poly> berlekamp_factor(const Poly& f) {
    int m = poly_degree(f);
    if (m <= 1) return {f};

    vector<Poly> rows(m, Poly{0, 0});
    Poly x_poly{2, 0};
    Poly x2 = poly_mul_mod(x_poly, x_poly, f);
    Poly power{1, 0};

    // Build Q matrix columns from x^(2j) mod f, then compute Q - I.
    for (int col = 0; col < m; ++col) {
        for (int row = 0; row < m; ++row) {
            if (poly_get_bit(power, row)) {
                poly_toggle_bit(rows[row], col);
            }
        }
        power = poly_mul_mod(power, x2, f);
    }
    for (int i = 0; i < m; ++i) poly_toggle_bit(rows[i], i);

    // Row-reduce to RREF.
    vector<int> pivot_col;
    pivot_col.reserve(m);
    int rank = 0;
    for (int col = 0; col < m; ++col) {
        int sel = -1;
        for (int r = rank; r < m; ++r) {
            if (poly_get_bit(rows[r], col)) {
                sel = r;
                break;
            }
        }
        if (sel == -1) continue;
        swap(rows[rank], rows[sel]);
        pivot_col.push_back(col);
        for (int r = 0; r < m; ++r) {
            if (r != rank && poly_get_bit(rows[r], col)) {
                rows[r] ^= rows[rank];
            }
        }
        ++rank;
    }

    vector<bool> is_pivot(m, false);
    for (int i = 0; i < rank; ++i) is_pivot[pivot_col[i]] = true;

    vector<Poly> basis;
    basis.reserve(m - rank);
    for (int col = 0; col < m; ++col) {
        if (is_pivot[col]) continue;
        Poly vec{0, 0};
        poly_toggle_bit(vec, col);
        for (int i = 0; i < rank; ++i) {
            int p = pivot_col[i];
            if (poly_get_bit(rows[i], col)) poly_toggle_bit(vec, p);
        }
        basis.push_back(vec);
    }

    vector<Poly> factors;
    factors.push_back(f);
    for (const Poly& b : basis) {
        if (b.is_zero() || b.is_one()) continue;
        vector<Poly> next;
        next.reserve(factors.size() * 2);
        for (const Poly& h : factors) {
            if (poly_degree(h) <= 1) {
                next.push_back(h);
                continue;
            }
            Poly g = poly_gcd(h, b);
            if (!g.is_one() && !poly_equal(g, h)) {
                Poly q = poly_div(h, g);
                next.push_back(g);
                next.push_back(q);
                continue;
            }
            Poly b1 = b;
            poly_toggle_bit(b1, 0);
            g = poly_gcd(h, b1);
            if (!g.is_one() && !poly_equal(g, h)) {
                Poly q = poly_div(h, g);
                next.push_back(g);
                next.push_back(q);
            } else {
                next.push_back(h);
            }
        }
        factors.swap(next);
    }
    return factors;
}

uint64_t order_in_factor(const Poly& lambda, const Poly& mod, const vector<uint32_t>& primes) {
    if (lambda.is_zero()) return 0;
    if (lambda.is_one()) return 1;

    int deg = poly_degree(mod);
    Poly cur = lambda;
    int e = 0;
    do {
        cur = poly_mul_mod(cur, cur, mod);
        ++e;
    } while (!poly_equal(cur, lambda) && e <= deg);

    if (e > deg || e > 63) {
        cerr << "Validation failure: unexpected field degree " << e << '\n';
        exit(1);
    }

    uint64_t order = (e == 64) ? ~0ULL : ((1ULL << e) - 1ULL);
    uint64_t reduced = order;
    vector<uint64_t> factors = factorize_uint64(order, primes);
    for (uint64_t p : factors) {
        while (reduced % p == 0) {
            uint64_t trial = reduced / p;
            if (poly_pow_mod(lambda, trial, mod).is_one()) {
                reduced = trial;
            } else {
                break;
            }
        }
    }
    return reduced;
}

vector<uint64_t> compute_orders_for_m(int m, const vector<uint32_t>& primes) {
    Poly f{1, 0};
    poly_toggle_bit(f, m);

    vector<Poly> factors = berlekamp_factor(f);
    vector<uint64_t> orders;
    orders.reserve(factors.size());
    Poly x_poly{2, 0};
    for (const Poly& g : factors) {
        Poly x_pow = poly_pow_mod(x_poly, static_cast<uint64_t>(m - 1), g);
        Poly lambda = x_poly ^ x_pow;
        lambda = poly_mod(lambda, g);
        orders.push_back(order_in_factor(lambda, g, primes));
    }
    return orders;
}

uint64_t lcm_u64(uint64_t a, uint64_t b) {
    return a / std::gcd(a, b) * b;
}

unordered_set<uint64_t> periods_for_n(int n,
                                      const vector<uint64_t>& orders,
                                      int k) {
    unordered_set<uint64_t> periods;
    periods.insert(1);
    for (uint64_t ord : orders) {
        vector<uint64_t> opts;
        opts.reserve(static_cast<size_t>(k + 2));
        opts.push_back(1);
        if (ord == 0) {
            // lambda = 0 contributes only period 1.
        } else if (ord == 1) {
            for (int s = 1; s <= k; ++s) opts.push_back(1ULL << s);
        } else {
            for (int s = 0; s <= k; ++s) opts.push_back(ord << s);
        }

        unordered_set<uint64_t> next;
        next.reserve(periods.size() * opts.size());
        for (uint64_t p : periods) {
            for (uint64_t o : opts) {
                next.insert(lcm_u64(p, o));
            }
        }
        periods.swap(next);
    }
    return periods;
}

uint64_t compute_S(int N,
                   unordered_map<int, vector<uint64_t>>& order_cache,
                   const vector<uint32_t>& primes) {
    unordered_set<uint64_t> global_periods;
    for (int n = 3; n <= N; ++n) {
        int m = n;
        int k = 0;
        while ((m & 1) == 0) {
            m >>= 1;
            ++k;
        }
        auto it = order_cache.find(m);
        if (it == order_cache.end()) {
            auto orders = compute_orders_for_m(m, primes);
            it = order_cache.emplace(m, std::move(orders)).first;
        }
        unordered_set<uint64_t> per = periods_for_n(n, it->second, k);
        global_periods.insert(per.begin(), per.end());
    }
    uint64_t sum = 0;
    for (uint64_t v : global_periods) sum += v;
    return sum;
}

bool check_periods(int n,
                   const vector<uint64_t>& expected,
                   unordered_map<int, vector<uint64_t>>& order_cache,
                   const vector<uint32_t>& primes) {
    int m = n;
    int k = 0;
    while ((m & 1) == 0) {
        m >>= 1;
        ++k;
    }
    auto it = order_cache.find(m);
    if (it == order_cache.end()) {
        auto orders = compute_orders_for_m(m, primes);
        it = order_cache.emplace(m, std::move(orders)).first;
    }
    unordered_set<uint64_t> per = periods_for_n(n, it->second, k);
    vector<uint64_t> got(per.begin(), per.end());
    vector<uint64_t> exp = expected;
    sort(got.begin(), got.end());
    sort(exp.begin(), exp.end());
    return got == exp;
}

}  // namespace

int main() {
    const int N = 100;
    vector<uint32_t> primes = sieve_primes(5'000'000);

    unordered_map<int, vector<uint64_t>> order_cache;

    if (!check_periods(5, {1, 3}, order_cache, primes) ||
        !check_periods(6, {1, 2}, order_cache, primes)) {
        cerr << "Validation failure: period sets for n=5 or n=6 mismatch\n";
        return 1;
    }

    uint64_t s30 = compute_S(30, order_cache, primes);
    if (s30 != 20381ULL) {
        cerr << "Validation failure: S(30) mismatch\n";
        return 1;
    }

    uint64_t answer = compute_S(N, order_cache, primes);
    cout << answer << '\n';
    return 0;
}
