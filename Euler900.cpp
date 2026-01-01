#include <algorithm>
#include <cassert>
#include <cstdint>
#include <future>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

static constexpr long long MOD = 900497239LL;

static long long mod_mul(long long a, long long b) {
    return static_cast<long long>((__int128)a * b % MOD);
}

static long long mod_pow(long long base, long long exp) {
    long long result = 1 % MOD;
    base %= MOD;
    while (exp > 0) {
        if (exp & 1LL) {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        exp >>= 1LL;
    }
    return result;
}

static long long mod_inv(long long a) {
    long long b = MOD;
    long long x0 = 1, y0 = 0;
    long long x1 = 0, y1 = 1;
    while (b != 0) {
        long long q = a / b;
        long long t = a % b;
        a = b;
        b = t;
        t = x0 - q * x1;
        x0 = x1;
        x1 = t;
        t = y0 - q * y1;
        y0 = y1;
        y1 = t;
    }
    if (a != 1) {
        cerr << "Modular inverse does not exist.\n";
        exit(1);
    }
    long long res = x0 % MOD;
    if (res < 0) res += MOD;
    return res;
}

static long long compute_S_closed(long long N) {
    long long inv3 = mod_inv(3);
    long long inv7 = mod_inv(7);

    long long pow4N = mod_pow(4, N);
    long long sumA = mod_mul((pow4N - 1 + MOD) % MOD, inv3);

    long long sumC = (mod_pow(2, N + 1) - 2) % MOD;
    if (sumC < 0) sumC += MOD;

    long long K = N / 2;
    long long pow8K = mod_pow(8, K);
    long long geo = mod_mul((pow8K - 1 + MOD) % MOD, inv7);
    long long sumB = mod_mul(3, geo);
    if (N % 2 == 1) {
        sumB += pow8K;
        if (sumB >= MOD) sumB -= MOD;
    }

    long long ans = (sumA + sumB - sumC) % MOD;
    if (ans < 0) ans += MOD;
    return ans;
}

static inline int bitlen_u64(uint64_t x) {
    return 64 - __builtin_clzll(x);
}

static uint64_t t_formula_u64(uint64_t n) {
    int m = bitlen_u64(n);
    uint64_t mask = (m == 64) ? ~0ULL : ((1ULL << m) - 1ULL);
    uint64_t r = static_cast<uint64_t>((__int128)n * n) & mask;
    uint64_t target = (n & 1ULL) ? mask : 0ULL;
    return (target - r) & mask;
}

static long long compute_S_direct_small(int N) {
    uint64_t lim = 1ULL << N;
    __int128 acc = 0;
    for (uint64_t n = 1; n <= lim; ++n) {
        acc += t_formula_u64(n);
    }
    return static_cast<long long>(acc % MOD);
}

struct VecHash {
    size_t operator()(const vector<int> &v) const noexcept {
        size_t h = 1469598103934665603ULL;
        for (int x : v) {
            h ^= static_cast<size_t>(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

struct BruteSolver {
    unordered_map<vector<int>, bool, VecHash> memo;

    bool is_losing(vector<int> piles) {
        sort(piles.begin(), piles.end());
        auto it = memo.find(piles);
        if (it != memo.end()) return it->second;

        int m = static_cast<int>(piles.size());
        int s = piles[0];

        long long maxSum = 0;
        for (int a : piles) {
            maxSum += (a - 1);
        }
        if (maxSum < s) {
            memo[piles] = true;
            return true;
        }

        vector<int> take(m, 0);
        function<bool(int, int)> dfs = [&](int i, int remaining) -> bool {
            if (i == m) {
                if (remaining != 0) return false;
                vector<int> next = piles;
                for (int j = 0; j < m; ++j) {
                    next[j] -= take[j];
                }
                if (is_losing(next)) return true;
                return false;
            }
            int ai = piles[i];
            int up = min(ai - 1, remaining);
            for (int x = 0; x <= up; ++x) {
                take[i] = x;
                if (dfs(i + 1, remaining - x)) return true;
            }
            take[i] = 0;
            return false;
        };

        bool hasWinningMove = dfs(0, s);
        bool res = !hasWinningMove;
        memo[piles] = res;
        return res;
    }

    int t_bruteforce(int n) {
        int m = 32 - __builtin_clz(n);
        int M = 1 << m;
        for (int k = 0; k < M; ++k) {
            vector<int> piles(n + 1, n);
            piles.back() = n + k;
            if (is_losing(piles)) return k;
        }
        return -1;
    }
};

static void run_validations() {
    {
        long long got = compute_S_closed(10);
        if (got != 361522) {
            cerr << "Validation failed: S(10) expected 361522, got " << got << "\n";
            exit(1);
        }
    }

    for (int N = 1; N <= 14; ++N) {
        long long a = compute_S_closed(N);
        long long b = compute_S_direct_small(N);
        if (a != b) {
            cerr << "Validation failed: N=" << N << " closed=" << a << " direct=" << b << "\n";
            exit(1);
        }
    }

    vector<future<void>> futs;
    for (int n = 1; n <= 8; ++n) {
        futs.push_back(async(launch::async, [n]() {
            BruteSolver bs;
            int brute = bs.t_bruteforce(n);
            uint64_t pred = t_formula_u64(static_cast<uint64_t>(n));
            if (static_cast<uint64_t>(brute) != pred) {
                cerr << "Validation failed: t(" << n << ") brute=" << brute
                     << " pred=" << pred << "\n";
                exit(1);
            }
        }));
    }
    for (auto &f : futs) {
        f.get();
    }
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    bool validate = false;
    long long N = 10000;
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        if (s == "--validate") {
            validate = true;
        } else {
            N = stoll(s);
        }
    }

    if (validate) {
        run_validations();
    }

    cout << compute_S_closed(N) << "\n";
    return 0;
}
