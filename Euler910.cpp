#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>

using namespace std;

static constexpr int DEG = 40;
static constexpr uint32_t MOD = 1000000000u;

struct Poly {
    array<uint32_t, DEG> c{};

    uint32_t eval(uint32_t x) const {
        uint64_t res = 0;
        uint64_t p = 1;
        for (int i = 0; i < DEG; ++i) {
            res += (p * c[i]) % MOD;
            res %= MOD;
            p = (p * x) % MOD;
        }
        return static_cast<uint32_t>(res);
    }
};

static Poly exponent_poly;

static inline Poly add_poly(Poly a, const Poly& b) {
    for (int i = 0; i < DEG; ++i) {
        uint32_t v = a.c[i] + b.c[i];
        if (v >= MOD) {
            v -= MOD;
        }
        a.c[i] = v;
    }
    return a;
}

static inline Poly lambda_poly(Poly a, uint32_t k) {
    for (int i = 0; i < DEG; ++i) {
        a.c[i] = static_cast<uint32_t>((static_cast<uint64_t>(a.c[i]) * k) % MOD);
    }
    return a;
}

static inline Poly shift_poly(Poly p) {
    uint32_t last = p.c[DEG - 1];
    for (int i = DEG - 1; i > 0; --i) {
        p.c[i] = p.c[i - 1];
    }
    p.c[0] = 0;
    return add_poly(p, lambda_poly(exponent_poly, last));
}

static inline Poly mul_poly(const Poly& p, Poly q) {
    Poly res;
    for (int i = 0; i < DEG; ++i) {
        res = add_poly(res, lambda_poly(q, p.c[i]));
        q = shift_poly(q);
    }
    return res;
}

static Poly pow_poly(Poly p, uint32_t n) {
    Poly res;
    res.c[0] = 1;
    while (n > 0) {
        if (n & 1u) {
            res = mul_poly(res, p);
        }
        n >>= 1u;
        if (n > 0) {
            p = mul_poly(p, p);
        }
    }
    return res;
}

static inline Poly compose_poly(const Poly& p, const Poly& q) {
    Poly res;
    Poly prod;
    prod.c[0] = 1;
    for (int i = 0; i < DEG; ++i) {
        res = add_poly(res, lambda_poly(prod, p.c[i]));
        prod = mul_poly(prod, q);
    }
    return res;
}

static Poly iter_compose(Poly p, uint32_t n) {
    if (n == 0) {
        Poly id;
        id.c[0] = 1;
        return id;
    }

    Poly res;
    bool init = false;
    while (n > 0) {
        if (n & 1u) {
            if (!init) {
                res = p;
                init = true;
            } else {
                res = compose_poly(p, res);
            }
        }
        n >>= 1u;
        if (n > 0) {
            p = compose_poly(p, p);
        }
    }
    return res;
}

static inline Poly d1(uint32_t n) {
    Poly p;
    p.c[1] = 1;
    Poly q = pow_poly(p, n);
    return add_poly(q, shift_poly(q));
}

static inline Poly d2(uint32_t n, const Poly& u) {
    Poly w = iter_compose(u, n);
    return compose_poly(w, shift_poly(u));
}

static Poly d3(uint32_t n, uint32_t m, uint32_t l) {
    if (n == 0) {
        Poly u = d1(l);
        Poly v = d2(m, u);
        return compose_poly(u, v);
    }
    Poly u = d3(n - 1, m, l);
    return d2(m, u);
}

static void init_exponent_poly() {
    exponent_poly.c.fill(0);
    exponent_poly.c[0] = 1;

    for (uint32_t i = 1; i <= DEG; ++i) {
        for (int j = DEG - 1; j > 0; --j) {
            exponent_poly.c[j] = static_cast<uint32_t>((exponent_poly.c[j - 1] + static_cast<uint64_t>(i) * exponent_poly.c[j]) % MOD);
        }
        exponent_poly.c[0] = static_cast<uint32_t>((static_cast<uint64_t>(exponent_poly.c[0]) * i) % MOD);
    }

    for (int i = 0; i < DEG; ++i) {
        exponent_poly.c[i] = (MOD - exponent_poly.c[i]) % MOD;
    }
}

static uint32_t solve(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e) {
    Poly p = d3(a, b, c);
    return (p.eval(d) + e) % MOD;
}

static bool run_validation() {
    struct T {
        uint32_t a;
        uint32_t b;
        uint32_t c;
        uint32_t d;
        uint32_t expected;
    };

    const T tests[] = {
        {0, 1, 1, 1, 42},
        {0, 2, 1, 2, 599882556},
        {1, 1, 1, 2, 707063360},
        {1, 2, 2, 3, 600000000},
        {2, 1, 2, 2, 0},
    };

    for (const auto& t : tests) {
        uint32_t got = solve(t.a, t.b, t.c, t.d, 0);
        if (got != t.expected) {
            cerr << "Validation failed for a=" << t.a << " b=" << t.b << " c=" << t.c
                 << " d=" << t.d << " expected=" << t.expected << " got=" << got << '\n';
            return false;
        }
    }

    return true;
}

int main() {
    const uint32_t a = 12;
    const uint32_t b = 345678;
    const uint32_t c = 9012345;
    const uint32_t d = 678;
    const uint32_t e = 90;

    init_exponent_poly();

    if (!run_validation()) {
        return 1;
    }

    uint32_t ans = solve(a, b, c, d, e);
    cout << setw(9) << setfill('0') << ans << '\n';
    return 0;
}
