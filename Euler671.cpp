#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

using ll = long long;

static const ll MOD = 1000004321LL;

struct State {
    int rt;
    int rb;
    int pv;
    int ct;
    int cb;
};

struct MatrixData {
    int k;
    int S;
    vector<State> states;
    vector<int> id;
    vector<vector<int>> adj;
};

struct RecurrenceData {
    vector<ll> coef;
    vector<ll> init;
};

int encode_state(int rt, int rb, int pv, int ct, int cb, int k) {
    return ((((rt * 3 + rb) * 2 + pv) * k + ct) * k + cb);
}

MatrixData build_matrix(int k) {
    MatrixData md;
    md.k = k;
    int total = 3 * 3 * 2 * k * k;
    md.id.assign(total, -1);
    md.states.reserve(total);

    for (int rt = 0; rt <= 2; ++rt) {
        for (int rb = 0; rb <= 2; ++rb) {
            for (int pv = 0; pv <= 1; ++pv) {
                for (int ct = 0; ct < k; ++ct) {
                    for (int cb = 0; cb < k; ++cb) {
                        bool ok = false;
                        if (pv == 1) {
                            ok = (rt == 0 && rb == 0 && ct == cb);
                        } else {
                            ok = (ct != cb);
                        }
                        if (!ok) continue;
                        int idx = static_cast<int>(md.states.size());
                        md.states.push_back({rt, rb, pv, ct, cb});
                        md.id[encode_state(rt, rb, pv, ct, cb, k)] = idx;
                    }
                }
            }
        }
    }

    md.S = static_cast<int>(md.states.size());
    md.adj.assign(md.S, {});

    auto idx_of = [&](int rt, int rb, int pv, int ct, int cb) {
        return md.id[encode_state(rt, rb, pv, ct, cb, k)];
    };

    for (int i = 0; i < md.S; ++i) {
        const State& s = md.states[i];
        if (s.rt > 0 && s.rb > 0) {
            md.adj[i].push_back(idx_of(s.rt - 1, s.rb - 1, 0, s.ct, s.cb));
            continue;
        }
        if (s.rt > 0 && s.rb == 0) {
            for (int lb = 1; lb <= 3; ++lb) {
                for (int cb = 0; cb < k; ++cb) {
                    if (cb == s.ct || cb == s.cb) continue;
                    md.adj[i].push_back(idx_of(s.rt - 1, lb - 1, 0, s.ct, cb));
                }
            }
            continue;
        }
        if (s.rt == 0 && s.rb > 0) {
            for (int lt = 1; lt <= 3; ++lt) {
                for (int ct = 0; ct < k; ++ct) {
                    if (ct == s.ct || ct == s.cb) continue;
                    md.adj[i].push_back(idx_of(lt - 1, s.rb - 1, 0, ct, s.cb));
                }
            }
            continue;
        }

        for (int c = 0; c < k; ++c) {
            if (c == s.ct || c == s.cb) continue;
            md.adj[i].push_back(idx_of(0, 0, 1, c, c));
        }

        if (s.pv == 1) {
            // Avoid four-corners: separate starts only after a vertical domino.
            for (int lt = 1; lt <= 3; ++lt) {
                for (int lb = 1; lb <= 3; ++lb) {
                    for (int ct = 0; ct < k; ++ct) {
                        if (ct == s.ct) continue;
                        for (int cb = 0; cb < k; ++cb) {
                            if (cb == s.cb || cb == ct) continue;
                            md.adj[i].push_back(idx_of(lt - 1, lb - 1, 0, ct, cb));
                        }
                    }
                }
            }
        }
    }

    return md;
}

using Matrix = vector<vector<uint64_t>>;

int choose_threads(int S) {
    unsigned hc = thread::hardware_concurrency();
    if (hc <= 1 || S < 200) return 1;
    int threads = static_cast<int>(hc);
    threads = min(threads, 8);
    threads = min(threads, S);
    return max(threads, 1);
}

Matrix multiply_dense_sparse(const Matrix& P, const vector<vector<int>>& adj, int threads) {
    int S = static_cast<int>(P.size());
    Matrix next(S, vector<uint64_t>(S, 0));

    auto worker = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            const auto& row = P[i];
            auto& out = next[i];
            for (int k = 0; k < S; ++k) {
                uint64_t val = row[k];
                if (val == 0) continue;
                const auto& edges = adj[k];
                for (int to : edges) {
                    out[to] += val;
                }
            }
            for (int j = 0; j < S; ++j) {
                out[j] %= MOD;
            }
        }
    };

    if (threads <= 1) {
        worker(0, S);
        return next;
    }

    vector<thread> pool;
    int block = (S + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
        int start = t * block;
        int end = min(S, start + block);
        if (start >= end) break;
        pool.emplace_back(worker, start, end);
    }
    for (auto& th : pool) th.join();
    return next;
}

ll mod_pow(ll base, ll exp) {
    ll res = 1 % MOD;
    ll cur = base % MOD;
    while (exp > 0) {
        if (exp & 1) res = static_cast<ll>((__int128)res * cur % MOD);
        cur = static_cast<ll>((__int128)cur * cur % MOD);
        exp >>= 1;
    }
    return res;
}

ll mod_inv(ll x) {
    return mod_pow((x % MOD + MOD) % MOD, MOD - 2);
}

vector<ll> berlekamp_massey(const vector<ll>& s) {
    vector<ll> C(1, 1), B(1, 1);
    ll L = 0;
    ll m = 1;
    ll b = 1;
    for (int n = 0; n < static_cast<int>(s.size()); ++n) {
        ll d = 0;
        for (int i = 0; i <= L; ++i) {
            d = (d + C[i] * s[n - i]) % MOD;
        }
        if (d == 0) {
            ++m;
            continue;
        }
        vector<ll> T = C;
        ll coef = d * mod_inv(b) % MOD;
        if (C.size() < B.size() + m) C.resize(B.size() + m, 0);
        for (int i = 0; i < static_cast<int>(B.size()); ++i) {
            ll val = (C[i + m] + MOD - coef * B[i] % MOD) % MOD;
            C[i + m] = val;
        }
        if (2 * L <= n) {
            L = n + 1 - L;
            B = T;
            b = d;
            m = 1;
        } else {
            ++m;
        }
    }
    C.erase(C.begin());
    for (ll& x : C) x = (MOD - x) % MOD;
    return C;
}

vector<ll> combine_poly(const vector<ll>& a, const vector<ll>& b, const vector<ll>& coef) {
    int L = static_cast<int>(coef.size());
    vector<ll> res(2 * L, 0);
    for (int i = 0; i < L; ++i) {
        if (a[i] == 0) continue;
        for (int j = 0; j < L; ++j) {
            if (b[j] == 0) continue;
            res[i + j] = (res[i + j] + a[i] * b[j]) % MOD;
        }
    }
    for (int i = 2 * L - 2; i >= L; --i) {
        if (res[i] == 0) continue;
        ll val = res[i];
        for (int j = 1; j <= L; ++j) {
            res[i - j] = (res[i - j] + val * coef[j - 1]) % MOD;
        }
    }
    res.resize(L);
    return res;
}

ll linear_rec(const vector<ll>& init, const vector<ll>& coef, long long n) {
    int L = static_cast<int>(coef.size());
    if (L == 0) return 0;
    if (n < static_cast<long long>(init.size())) return init[n];
    vector<ll> pol(L, 0), e(L, 0);
    pol[0] = 1;
    if (L == 1) {
        e[0] = coef[0];
    } else {
        e[1] = 1;
    }
    long long m = n;
    while (m > 0) {
        if (m & 1) pol = combine_poly(pol, e, coef);
        e = combine_poly(e, e, coef);
        m >>= 1;
    }
    ll res = 0;
    for (int i = 0; i < L; ++i) {
        res = (res + pol[i] * init[i]) % MOD;
    }
    return res;
}

RecurrenceData build_recurrence(int k) {
    MatrixData md = build_matrix(k);
    int S = md.S;
    int threads = choose_threads(S);

    Matrix P(S, vector<uint64_t>(S, 0));
    for (int i = 0; i < S; ++i) P[i][i] = 1;

    vector<ll> seq;
    seq.reserve(512);
    seq.push_back(S % MOD);

    vector<ll> coef;
    int last_L = 0;
    int stable = 0;
    int max_steps = 2000;

    for (int step = 1; step <= max_steps; ++step) {
        P = multiply_dense_sparse(P, md.adj, threads);
        ll tr = 0;
        for (int i = 0; i < S; ++i) {
            tr += static_cast<ll>(P[i][i]);
            if (tr >= MOD) tr %= MOD;
        }
        seq.push_back(tr % MOD);

        coef = berlekamp_massey(seq);
        int L = static_cast<int>(coef.size());
        if (L == last_L) {
            ++stable;
        } else {
            stable = 0;
            last_L = L;
        }
        if (L > 0 && static_cast<int>(seq.size()) >= 2 * L + 5 && stable >= 5) {
            break;
        }
    }

    if (coef.empty()) {
        coef.push_back(0);
    }

    vector<ll> init(coef.size(), 0);
    for (size_t i = 0; i < init.size() && i < seq.size(); ++i) {
        init[i] = seq[i];
    }

    return {coef, init};
}

vector<pair<ll, int>> factorize(ll n) {
    vector<pair<ll, int>> factors;
    for (ll p = 2; p * p <= n; p += (p == 2 ? 1 : 2)) {
        if (n % p != 0) continue;
        int cnt = 0;
        while (n % p == 0) {
            n /= p;
            ++cnt;
        }
        factors.push_back({p, cnt});
    }
    if (n > 1) factors.push_back({n, 1});
    return factors;
}

void gen_divisors(int idx, ll cur_d, ll cur_phi,
                  const vector<pair<ll, int>>& factors,
                  vector<pair<ll, ll>>& out) {
    if (idx == static_cast<int>(factors.size())) {
        out.push_back({cur_d, cur_phi % MOD});
        return;
    }
    ll p = factors[idx].first;
    int e = factors[idx].second;
    gen_divisors(idx + 1, cur_d, cur_phi, factors, out);

    ll p_pow = 1;
    ll phi_pow = 1;
    for (int i = 1; i <= e; ++i) {
        p_pow *= p;
        if (i == 1) {
            phi_pow = p - 1;
        } else {
            phi_pow *= p;
        }
        gen_divisors(idx + 1, cur_d * p_pow, cur_phi * phi_pow, factors, out);
    }
}

ll solve_problem(int k, ll n) {
    RecurrenceData rec = build_recurrence(k);
    auto factors = factorize(n);
    vector<pair<ll, ll>> divisors;
    gen_divisors(0, 1, 1, factors, divisors);

    ll total = 0;
    for (const auto& dv : divisors) {
        ll d = dv.first;
        ll phi_mod = dv.second % MOD;
        ll term = linear_rec(rec.init, rec.coef, n / d);
        total = (total + phi_mod * term) % MOD;
    }

    ll inv_n = mod_inv(n % MOD);
    return total * inv_n % MOD;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (solve_problem(4, 3) != 104) {
        cerr << "Validation failed: F4(3) != 104\n";
        return 1;
    }
    if (solve_problem(5, 7) != 3327300) {
        cerr << "Validation failed: F5(7) != 3327300\n";
        return 1;
    }
    if (solve_problem(6, 101) != 75309980) {
        cerr << "Validation failed: F6(101) != 75309980\n";
        return 1;
    }

    const ll n = 10004003002001LL;
    cout << solve_problem(10, n) << "\n";
    return 0;
}
