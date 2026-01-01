#include <bits/stdc++.h>

using namespace std;

static const uint32_t MOD = 1000000007u;

static inline uint32_t addmod(uint32_t a, uint32_t b){ uint32_t r = a + b; if(r >= MOD) r -= MOD; return r; }
static inline uint32_t submod(uint32_t a, uint32_t b){ return a >= b ? a - b : a + MOD - b; }
static inline uint32_t fastmod(uint64_t x){
    static const uint64_t inv = (uint64_t)((__uint128_t(1) << 64) / MOD);
    uint64_t q = (uint64_t)((__uint128_t)x * inv >> 64);
    uint64_t r = x - q * MOD;
    if(r >= MOD) r -= MOD;
    return (uint32_t)r;
}
static inline uint32_t mulmod(uint32_t a, uint32_t b){ return fastmod((uint64_t)a * b); }

static uint32_t modpow(uint32_t a, uint32_t e){
    uint32_t r = 1u;
    a %= MOD;
    while(e > 0){
        if(e & 1) r = mulmod(r, a);
        a = mulmod(a, a);
        e >>= 1;
    }
    return r;
}

struct Comb {
    int maxI;
    vector<uint32_t> fact, invfact, pow2, inv;
    explicit Comb(int maxI_=0){ init(maxI_); }
    void init(int maxI_){
        maxI = max(0, maxI_);
        fact.assign(maxI+1, 1);
        invfact.assign(maxI+1, 1);
        pow2.assign(maxI+1, 1);
        inv.assign(maxI+1, 0);
        for(int i=1;i<=maxI;i++){
            fact[i] = mulmod(fact[i-1], (uint32_t)i);
            pow2[i] = addmod(pow2[i-1], pow2[i-1]);
        }
        invfact[maxI] = modpow(fact[maxI], MOD-2);
        for(int i=maxI;i>=1;i--){
            invfact[i-1] = mulmod(invfact[i], (uint32_t)i);
        }
        if(maxI >= 1){
            inv[1] = 1;
            for(int i=2;i<=maxI;i++){
                inv[i] = MOD - mulmod((uint32_t)(MOD / i), inv[MOD % i]);
            }
        }
    }
};

static inline uint32_t coefP(const Comb& C, int i, int u, int v, int w){
    int n = u + v + w;
    int k = n - 2*i;
    if(k < 0 || k > i) return 0;
    if(u < k || v < k || w < k) return 0;
    if(((u-k)&1) || ((v-k)&1) || ((w-k)&1)) return 0;
    int p = (u-k)/2, q = (v-k)/2, r = (w-k)/2;
    if(p < 0 || q < 0 || r < 0) return 0;
    uint32_t res = C.fact[i];
    res = mulmod(res, C.invfact[k]);
    res = mulmod(res, C.invfact[p]);
    res = mulmod(res, C.invfact[q]);
    res = mulmod(res, C.invfact[r]);
    res = mulmod(res, C.pow2[k]);
    return res;
}

static uint32_t H(const Comb& C, int u, int v, int w){
    if(u < 0 || v < 0 || w < 0) return 0;
    int n = u + v + w;
    if(((u ^ n) & 1) || ((v ^ n) & 1) || ((w ^ n) & 1)) return 0;
    int A = (v + w) / 2;
    int B = (u + w) / 2;
    int Cc = (u + v) / 2;
    int lo = (n + 2) / 3;
    if(A > lo) lo = A;
    if(B > lo) lo = B;
    if(Cc > lo) lo = Cc;
    int hi = n / 2;
    if(lo > hi) return 0;
    int i0 = lo;
    int k0 = n - 2 * i0;
    int p0 = i0 - A;
    int q0 = i0 - B;
    int r0 = i0 - Cc;
    uint32_t coef = C.fact[i0];
    coef = mulmod(coef, C.invfact[k0]);
    coef = mulmod(coef, C.invfact[p0]);
    coef = mulmod(coef, C.invfact[q0]);
    coef = mulmod(coef, C.invfact[r0]);
    coef = mulmod(coef, C.pow2[k0]);
    uint32_t s = coef;
    const uint32_t inv4 = C.inv[4];
    int k = k0;
    for(int i = i0; i < hi; i++){
        uint32_t num = mulmod((uint32_t)(i + 1), mulmod((uint32_t)k, (uint32_t)(k - 1)));
        uint32_t den = mulmod(C.inv[i + 1 - A], mulmod(C.inv[i + 1 - B], C.inv[i + 1 - Cc]));
        coef = mulmod(coef, num);
        coef = mulmod(coef, den);
        coef = mulmod(coef, inv4);
        s = addmod(s, coef);
        k -= 2;
    }
    return s;
}

static uint32_t gCoeff(const Comb& C, int a, int b, int c){
    uint32_t res = 0;
    res = addmod(res, H(C, a,   b,   c));
    res = addmod(res, H(C, a-1, b,   c));
    res = addmod(res, H(C, a,   b,   c-1));
    res = addmod(res, H(C, a-1, b-1, c));
    res = addmod(res, H(C, a,   b-1, c-1));
    res = submod(res, H(C, a,   b-2, c));
    return res;
}

static vector<array<int,3>> enumerateXorTriples(int n){
    vector<array<int,3>> cur;
    cur.push_back({0,0,0});
    if(n & 1) return {};
    for(int bit = 1; (1<<bit) <= n*2 || bit <= 20; bit++){
        if((n >> bit) & 1){
            int v = 1 << (bit-1);
            vector<array<int,3>> nxt;
            nxt.reserve(cur.size() * 3);
            for(auto t: cur){
                nxt.push_back({t[0], t[1] + v, t[2] + v});
                nxt.push_back({t[0] + v, t[1], t[2] + v});
                nxt.push_back({t[0] + v, t[1] + v, t[2]});
            }
            cur.swap(nxt);
        }
        if((1<<bit) > n && bit > 20) break;
    }
    vector<array<int,3>> out;
    out.reserve(cur.size());
    for(auto t: cur){
        if(t[0] + t[1] + t[2] == n && ((t[0] ^ t[1] ^ t[2]) == 0)) out.push_back(t);
    }
    return out;
}

static uint32_t solveFast(const Comb& C, int n){
    if(n < 0) return 0;
    if(n & 1) return 0;
    auto triples = enumerateXorTriples(n);

    unsigned hw = thread::hardware_concurrency();
    unsigned T = max(1u, hw);
    T = min<unsigned>(T, (unsigned)triples.size());
    if(T == 0) T = 1;

    vector<uint32_t> partial(T, 0);
    vector<thread> threads;
    threads.reserve(T);
    atomic<size_t> done(0);
    atomic<bool> running(true);
    const size_t total = triples.size();
    thread progress;
    if(total > 0){
        progress = thread([&](){
            const int width = 40;
            auto start = chrono::steady_clock::now();
            while(running.load(memory_order_relaxed)){
                size_t d = done.load(memory_order_relaxed);
                double ratio = total ? (double)d / (double)total : 1.0;
                if(ratio > 1.0) ratio = 1.0;
                int filled = (int)(ratio * width);
                auto now = chrono::steady_clock::now();
                double secs = chrono::duration<double>(now - start).count();
                cerr << "\r[";
                for(int i=0;i<width;i++) cerr << (i < filled ? '#' : '-');
                cerr << "] " << fixed << setprecision(2) << (ratio * 100.0) << "% "
                     << d << "/" << total << " " << (uint64_t)secs << "s" << flush;
                this_thread::sleep_for(chrono::milliseconds(200));
            }
        });
    }

    auto worker = [&](unsigned tid){
        size_t L = triples.size();
        size_t from = (L * tid) / T;
        size_t to   = (L * (tid + 1)) / T;
        uint32_t acc = 0;
        for(size_t i = from; i < to; i++){
            auto t = triples[i];
            acc = addmod(acc, gCoeff(C, t[0], t[1], t[2]));
            done.fetch_add(1, memory_order_relaxed);
        }
        partial[tid] = acc;
    };
    for(unsigned t=0; t<T; t++) threads.emplace_back(worker, t);
    for(auto& th: threads) th.join();
    running.store(false, memory_order_relaxed);
    if(progress.joinable()){
        progress.join();
        cerr << "\r";
        cerr << string(80, ' ') << "\r";
    }

    uint32_t K = 0;
    for(unsigned t=0; t<T; t++) K = addmod(K, partial[t]);

    const uint32_t inv2 = (MOD + 1) / 2;
    uint32_t pow2n = modpow(2u, (uint32_t)n);
    uint32_t ans = mulmod(K, submod(pow2n, 1));
    ans = mulmod(ans, inv2);
    return ans;
}

static uint32_t solveBrute(int n){
    if(n < 0) return 0;
    if(n == 0) return 0;
    int N = n;

    vector<uint8_t> loc(N+1, 0);
    int cnt0 = N, cnt1 = 0, cnt2 = 0;
    vector<int8_t> dir(N+1, 0);
    for(int d=1; d<=N; d++){
        dir[d] = (((N-d)&1) ? 1 : -1);
    }

    uint64_t sumIdx = 0;
    auto isLose = [&](int a,int b,int c){ return ((a ^ b ^ c) == 0); };
    if(isLose(cnt0,cnt1,cnt2)) sumIdx = 0;
    uint32_t totalStates = (N >= 31 ? 0u : (1u << N));

    for(uint32_t m = 1; m < totalStates; m++){
        int d = __builtin_ctz(m) + 1;
        int oldp = loc[d];
        int newp = oldp + dir[d];
        if(newp < 0) newp += 3;
        else if(newp >= 3) newp -= 3;
        loc[d] = (uint8_t)newp;
        if(oldp == 0) cnt0--;
        else if(oldp == 1) cnt1--;
        else cnt2--;
        if(newp == 0) cnt0++;
        else if(newp == 1) cnt1++;
        else cnt2++;

        if(isLose(cnt0,cnt1,cnt2)){
            sumIdx += (uint64_t)m;
            sumIdx %= MOD;
        }
    }
    return (uint32_t)sumIdx;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 100000;
    if(argc > 1){
        n = atoi(argv[1]);
    }

    int maxN = max(n, 10);
    int maxI = maxN / 2 + 5;
    Comb C(maxI);

    {
        uint32_t f4_fast = solveFast(C, 4);
        uint32_t f10_fast = solveFast(C, 10);
        if(f4_fast != 30 || f10_fast != 67518){
            cerr << "Validation failed (fast): f(4)=" << f4_fast
                 << ", f(10)=" << f10_fast << "\n";
            return 1;
        }
        uint32_t f4_brute = solveBrute(4);
        uint32_t f10_brute = solveBrute(10);
        if(f4_brute != 30 || f10_brute != 67518 || f4_brute != f4_fast || f10_brute != f10_fast){
            cerr << "Validation failed (brute): f(4)=" << f4_brute
                 << ", f(10)=" << f10_brute << "\n";
            return 1;
        }
    }

    cout << solveFast(C, n) << "\n";
    return 0;
}
