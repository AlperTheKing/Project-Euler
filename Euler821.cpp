#include <bits/stdc++.h>
using namespace std;

using int64 = long long;
using i128 = __int128_t;

static inline int64 count_coprime6(int64 x) {
    if (x <= 0) return 0;
    return x - x / 2 - x / 3 + x / 6;
}

static vector<int64> generate_smooth(int64 N) {
    vector<int64> v;
    for (int64 p2 = 1; p2 <= N; ) {
        for (int64 cur = p2; cur <= N; ) {
            v.push_back(cur);
            if (cur > N / 3) break;
            cur *= 3;
        }
        if (p2 > N / 2) break;
        p2 *= 2;
    }
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
    return v;
}

static vector<int64> generate_bad(int64 N) {
    vector<int64> b;
    auto push_if = [&](int64 x) { if (x <= N) b.push_back(x); };

    push_if(6);
    push_if(24);
    push_if(54);

    {
        int64 x = 384;
        while (x <= N) {
            b.push_back(x);
            if (x > N / 8) break;
            x *= 8;
        }
    }
    {
        int64 x = 243;
        while (x <= N) {
            b.push_back(x);
            if (x > N / 27) break;
            x *= 27;
        }
    }

    sort(b.begin(), b.end());
    b.erase(unique(b.begin(), b.end()), b.end());
    return b;
}

static int64 solve_fast(int64 N) {
    vector<int64> smooth = generate_smooth(N);
    if (smooth.empty() || smooth.front() != 1) {
        throw runtime_error("smooth list must start with 1");
    }

    vector<int64> bad = generate_bad(N);

    for (int64 x : bad) {
        if (!binary_search(smooth.begin(), smooth.end(), x)) {
            cerr << "Bad number not found in smooth list: " << x << "\n";
            throw runtime_error("bad number missing from smooth list");
        }
    }

    vector<int64> fAt(smooth.size(), 0);
    size_t ib = 0;
    int64 f = 0;
    for (size_t i = 0; i < smooth.size(); i++) {
        if (ib < bad.size() && smooth[i] == bad[ib]) {
            ++ib;
        } else {
            ++f;
        }
        fAt[i] = f;
        if (fAt[i] > (int64)(i + 1)) {
            throw runtime_error("fAt exceeded i+1 (logic error)");
        }
    }
    if (ib != bad.size()) {
        throw runtime_error("Not all bad numbers were consumed in scan");
    }

    int numThreads = 1;
#ifdef USE_THREADS
    numThreads = max(1u, thread::hardware_concurrency());
#endif

    vector<i128> partial(numThreads, 0);

    auto worker = [&](int tid) {
        size_t n = smooth.size();
        size_t lo = (n * (size_t)tid) / (size_t)numThreads;
        size_t hi = (n * (size_t)(tid + 1)) / (size_t)numThreads;

        i128 sum = 0;
        for (size_t i = lo; i < hi; i++) {
            int64 L = smooth[i];
            int64 R = (i + 1 < n ? smooth[i + 1] - 1 : N);
            if (R > N) R = N;

            int64 mLow  = N / (R + 1) + 1;
            int64 mHigh = N / L;
            if (mLow > mHigh) continue;

            int64 cnt = count_coprime6(mHigh) - count_coprime6(mLow - 1);
            sum += (i128)fAt[i] * (i128)cnt;
        }
        partial[tid] = sum;
    };

    if (numThreads == 1) {
        worker(0);
    } else {
        vector<thread> th;
        th.reserve(numThreads);
        for (int t = 0; t < numThreads; t++) th.emplace_back(worker, t);
        for (auto &tt : th) tt.join();
    }

    i128 ans = 0;
    for (int t = 0; t < numThreads; t++) ans += partial[t];

    return (int64)ans;
}

#ifdef VALIDATE
static int64 brute_component_f_small(int64 L) {
    vector<pair<int,int>> pts;
    vector<int64> vals;
    for (int a = 0, p2 = 1; (int64)p2 <= L; a++) {
        for (int b = 0, cur = p2; (int64)cur <= L; b++) {
            pts.push_back({a,b});
            vals.push_back(cur);
            if (cur > L/3) break;
            cur *= 3;
        }
        if (p2 > L/2) break;
        p2 *= 2;
    }
    int n = (int)vals.size();
    if (n > 64) {
        throw runtime_error("Validation brute_component_f_small exceeded 64 nodes");
    }

    vector<int64> w(n);
    for (int i = 0; i < n; i++) {
        int64 v = vals[i];
        w[i] = 1 + (2*v <= L) + (3*v <= L);
    }

    unordered_map<long long, long long> idx;
    idx.reserve(n*2);
    for (int i = 0; i < n; i++) {
        long long key = ((long long)pts[i].first << 32) ^ (unsigned)pts[i].second;
        idx[key] = i;
    }

    vector<uint64_t> adj(n, 0);
    auto get = [&](int a,int b)->int {
        long long key = ((long long)a << 32) ^ (unsigned)b;
        auto it = idx.find(key);
        return (it == idx.end() ? -1 : (int)it->second);
    };

    for (int i = 0; i < n; i++) {
        int a = pts[i].first, b = pts[i].second;
        const int da[6] = {1,-1,0,0,1,-1};
        const int db[6] = {0,0,1,-1,-1,1};
        for (int k = 0; k < 6; k++) {
            int j = get(a + da[k], b + db[k]);
            if (j >= 0) adj[i] |= (1ULL << j);
        }
    }

    unordered_map<uint64_t, int64> memo;
    memo.reserve(1<<16);

    function<int64(uint64_t)> dp = [&](uint64_t mask)->int64 {
        if (!mask) return 0;
        auto it = memo.find(mask);
        if (it != memo.end()) return it->second;

        int i = __builtin_ctzll(mask);
        uint64_t m0 = mask & ~(1ULL << i);
        int64 best = dp(m0);
        int64 take = w[i] + dp(m0 & ~adj[i]);
        if (take > best) best = take;
        memo[mask] = best;
        return best;
    };

    uint64_t full = (n == 64 ? ~0ULL : ((1ULL << n) - 1ULL));
    return dp(full);
}

static int64 brute_F_small(int64 N) {
    unordered_map<int64, int64> memo_f;
    memo_f.reserve(N);
    auto f = [&](int64 L)->int64 {
        auto it = memo_f.find(L);
        if (it != memo_f.end()) return it->second;
        int64 val = brute_component_f_small(L);
        memo_f[L] = val;
        return val;
    };

    int64 ans = 0;
    for (int64 m = 1; m <= N; m++) {
        if (std::gcd<int64>(m, 6) != 1) continue;
        ans += f(N / m);
    }
    return ans;
}

static void run_validation() {
    {
        int64 N = 6;
        int64 fast = solve_fast(N);
        int64 brute = brute_F_small(N);
        assert(fast == brute);
        assert(fast == 5);
    }
    {
        int64 N = 20;
        int64 fast = solve_fast(N);
        int64 brute = brute_F_small(N);
        assert(fast == brute);
        assert(fast == 19);
    }

    vector<int64> tests = {1,2,3,4,5,7,8,9,10,50,100,200,400,1000,2000,5000};
    for (auto N : tests) {
        int64 fast = solve_fast(N);
        int64 brute = brute_F_small(N);
        if (fast != brute) {
            cerr << "Validation failed at N=" << N
                 << " fast=" << fast << " brute=" << brute << "\n";
            assert(false);
        }
    }
}
#endif

int main() {
#ifdef VALIDATE
    run_validation();
#endif

    const int64 N = 10000000000000000LL;
    cout << solve_fast(N) << "\n";
    return 0;
}
