#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace std;

using u128 = unsigned __int128;
using i128 = __int128_t;

static inline uint64_t mul_mod_u64(uint64_t a, uint64_t b, uint64_t mod) {
    return static_cast<uint64_t>((u128)a * (u128)b % (u128)mod);
}

static inline uint64_t sub_mod_u64(uint64_t a, uint64_t b, uint64_t mod) {
    return (a >= b) ? (a - b) : static_cast<uint64_t>((u128)mod + a - b);
}

static uint64_t modinv_u64(uint64_t a, uint64_t mod) {
    a %= mod;
    if (a == 0) throw runtime_error("modinv: a=0");

    i128 t = 0, newt = 1;
    i128 r = static_cast<i128>(mod), newr = static_cast<i128>(a);

    while (newr != 0) {
        i128 q = r / newr;
        i128 tmp_t = t - q * newt;
        t = newt;
        newt = tmp_t;
        i128 tmp_r = r - q * newr;
        r = newr;
        newr = tmp_r;
    }

    if (r != 1) throw runtime_error("modinv: not invertible");
    if (t < 0) t += static_cast<i128>(mod);
    return static_cast<uint64_t>(t);
}

static uint64_t pow_mod_u64(uint64_t a, uint64_t e, uint64_t mod) {
    uint64_t res = 1 % mod;
    a %= mod;
    while (e) {
        if (e & 1) res = mul_mod_u64(res, a, mod);
        a = mul_mod_u64(a, a, mod);
        e >>= 1;
    }
    return res;
}

static bool is_prime_trial(uint64_t p) {
    if (p < 2) return false;
    if (p % 2 == 0) return p == 2;
    for (uint64_t d = 3; d * d <= p; d += 2) {
        if (p % d == 0) return false;
    }
    return true;
}

static vector<uint64_t> pow_table(uint64_t p, int S) {
    vector<uint64_t> pp(S + 1);
    pp[0] = 1;
    for (int i = 1; i <= S; i++) pp[i] = pp[i - 1] * p;
    return pp;
}

static vector<uint64_t> compute_coeff_C(uint64_t p, int S) {
    const vector<uint64_t> pPow = pow_table(p, S);
    const uint64_t mod = pPow[S];
    const int N = (S - 1) / 2;
    if (N <= 0) return {};

    auto gcd_u64 = [](uint64_t a, uint64_t b) {
        while (b) {
            uint64_t t = a % b;
            a = b;
            b = t;
        }
        return a;
    };

    vector<uint64_t> sample_n;
    sample_n.reserve(N);
    for (uint64_t n = 1; sample_n.size() < static_cast<size_t>(N); ++n) {
        if (n % p != 0) sample_n.push_back(n);
    }

    vector<uint64_t> b(N);
    for (int i = 0; i < N; i++) {
        uint64_t n = sample_n[static_cast<size_t>(i)];
        uint64_t sum = 0;
        const uint64_t up = p * n;
        for (uint64_t m = 1; m <= up; m++) {
            if (m % p == 0) continue;
            sum += modinv_u64(m, mod);
            if (sum >= mod) sum %= mod;
        }
        b[n - 1] = sum % mod;
    }

    vector<vector<uint64_t>> M(N, vector<uint64_t>(N + 1, 0));
    for (int i = 0; i < N; i++) {
        uint64_t n = sample_n[static_cast<size_t>(i)];
        for (int k = 0; k < N; k++) {
            M[i][k] = pow_mod_u64(n, static_cast<uint64_t>(2 * (k + 1)), mod);
        }
        M[i][N] = b[i];
    }

    for (int col = 0; col < N; col++) {
        int pivot = -1;
        for (int r = col; r < N; r++) {
            if (gcd_u64(M[r][col], mod) == 1) {
                pivot = r;
                break;
            }
        }
        if (pivot == -1) {
            throw runtime_error("Coefficient solve failed: no invertible pivot.");
        }
        if (pivot != col) swap(M[pivot], M[col]);

        uint64_t invPivot = modinv_u64(M[col][col], mod);
        for (int c = col; c <= N; c++) {
            M[col][c] = mul_mod_u64(M[col][c], invPivot, mod);
        }

        for (int r = 0; r < N; r++) {
            if (r == col) continue;
            uint64_t factor = M[r][col];
            if (factor == 0) continue;
            for (int c = col; c <= N; c++) {
                uint64_t sub = mul_mod_u64(factor, M[col][c], mod);
                M[r][c] = sub_mod_u64(M[r][c], sub, mod);
            }
        }
    }

    vector<uint64_t> C(N);
    for (int i = 0; i < N; i++) C[i] = M[i][N] % mod;

    {
        uint64_t testn = sample_n.back() + 1;
        while (testn % p == 0) ++testn;
        uint64_t btest = 0;
        const uint64_t up = p * testn;
        for (uint64_t m = 1; m <= up; m++) {
            if (m % p == 0) continue;
            btest += modinv_u64(m, mod);
            if (btest >= mod) btest %= mod;
        }
        btest %= mod;

        uint64_t poly = 0;
        for (int k = 0; k < N; k++) {
            uint64_t term = mul_mod_u64(C[k],
                                        pow_mod_u64(testn, static_cast<uint64_t>(2 * (k + 1)), mod),
                                        mod);
            poly += term;
            if (poly >= mod) poly %= mod;
        }
        poly %= mod;

        if (btest != poly) {
            throw runtime_error("Coefficient validation failed: b_n mismatch.");
        }
    }

    return C;
}

struct Node {
    uint64_t n;
    uint64_t H;
};

static uint64_t max_Jp(uint64_t p, int S, bool use_threads) {
    if (p < 3 || (p % 2 == 0)) throw runtime_error("Expected an odd prime p>=3.");
    if (!is_prime_trial(p)) throw runtime_error("p is not prime.");

    const vector<uint64_t> pPow = pow_table(p, S);
    const uint64_t modS = pPow[S];
    const vector<uint64_t> C = compute_coeff_C(p, S);

    vector<uint64_t> Hsmall(p);
    Hsmall[0] = 0;
    for (uint64_t n = 1; n < p; n++) {
        Hsmall[n] = (Hsmall[n - 1] + modinv_u64(n, modS)) % modS;
    }

    vector<Node> cur;
    cur.reserve(p);
    uint64_t maxn = 0;
    for (uint64_t n = 1; n < p; n++) {
        if (Hsmall[n] % p == 0) {
            cur.push_back({n, Hsmall[n]});
            maxn = max(maxn, n);
        }
    }

    auto expand_one = [&](const Node& parent, uint64_t mod_r, uint64_t mod_next) -> vector<Node> {
        const uint64_t n = parent.n;
        const uint64_t hn = parent.H % mod_r;

        const uint64_t hn_div_p = (hn / p) % mod_next;

        uint64_t poly = 0;
        uint64_t nmod = n % mod_next;
        uint64_t n2 = mul_mod_u64(nmod, nmod, mod_next);
        uint64_t pow_n2k = n2;
        for (size_t k = 0; k < C.size(); k++) {
            uint64_t ck = C[k] % mod_next;
            if (ck) {
                poly += mul_mod_u64(ck, pow_n2k, mod_next);
                if (poly >= mod_next) poly %= mod_next;
            }
            pow_n2k = mul_mod_u64(pow_n2k, n2, mod_next);
        }
        poly %= mod_next;

        uint64_t hp_n = (hn_div_p + poly) % mod_next;

        vector<Node> children;
        children.reserve(p);

        if (hp_n % p == 0) children.push_back({p * n, hp_n});

        uint64_t h = hp_n;
        for (uint64_t k = 1; k < p; k++) {
            uint64_t denom = p * n + k;
            uint64_t invd = modinv_u64(denom, mod_next);
            h += invd;
            if (h >= mod_next) h %= mod_next;
            if (h % p == 0) children.push_back({p * n + k, h});
        }

        return children;
    };

    int r = S;
    while (!cur.empty() && r >= 2) {
        const uint64_t mod_r = pPow[r];
        const uint64_t mod_next = pPow[r - 1];

        vector<Node> nxt;

        if (!use_threads || cur.size() < 4) {
            for (const auto& node : cur) {
                vector<Node> kids = expand_one(node, mod_r, mod_next);
                for (auto& ch : kids) {
                    maxn = max(maxn, ch.n);
                    nxt.push_back(ch);
                }
            }
        } else {
            unsigned hw = thread::hardware_concurrency();
            if (hw == 0) hw = 2;
            unsigned T = min<unsigned>(hw, static_cast<unsigned>(cur.size()));

            vector<vector<Node>> parts(T);
            vector<uint64_t> localMax(T, 0);
            vector<thread> threads;
            threads.reserve(T);

            for (unsigned t = 0; t < T; t++) {
                size_t L = static_cast<size_t>(t * cur.size() / T);
                size_t R = static_cast<size_t>((t + 1) * cur.size() / T);
                threads.emplace_back([&, t, L, R]() {
                    uint64_t mx = 0;
                    vector<Node> out;
                    for (size_t i = L; i < R; i++) {
                        vector<Node> kids = expand_one(cur[i], mod_r, mod_next);
                        for (auto& ch : kids) {
                            mx = max(mx, ch.n);
                            out.push_back(ch);
                        }
                    }
                    localMax[t] = mx;
                    parts[t].swap(out);
                });
            }
            for (auto& th : threads) th.join();

            for (unsigned t = 0; t < T; t++) {
                maxn = max(maxn, localMax[t]);
                nxt.insert(nxt.end(), parts[t].begin(), parts[t].end());
            }
        }

        cur.swap(nxt);
        r--;
    }

    if (!cur.empty() && r < 2) {
        throw runtime_error("Precision S too low: increase S.");
    }

    return maxn;
}

static uint64_t compute_M(uint64_t p, int S, bool use_threads) {
    uint64_t maxJ = max_Jp(p, S, use_threads);
    return p * maxJ + (p - 1);
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    uint64_t p = 137;
    if (argc >= 2) {
        p = stoull(argv[1]);
    }

    const int S = 8;

    {
        uint64_t M7 = compute_M(7, S, false);
        if (M7 != 719102ULL) {
            cerr << "Validation failed: M(7) expected 719102, got " << M7 << "\n";
            return 1;
        }
    }

    uint64_t ans = compute_M(p, S, true);
    cout << ans << "\n";
    return 0;
}
