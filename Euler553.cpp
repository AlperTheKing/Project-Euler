#include <algorithm>
#include <cstdint>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {

using int64 = long long;

constexpr int kMod = 1'000'000'007;
constexpr int kModMinus1 = kMod - 1;

int64 mod_pow(int64 base, int64 exp, int64 mod) {
    int64 res = 1 % mod;
    int64 cur = base % mod;
    while (exp > 0) {
        if (exp & 1LL) res = static_cast<int64>((__int128)res * cur % mod);
        cur = static_cast<int64>((__int128)cur * cur % mod);
        exp >>= 1LL;
    }
    return res;
}

int mod_inv(int a) {
    return static_cast<int>(mod_pow(a, kMod - 2, kMod));
}

int mulmod(int64 a, int64 b) {
    return static_cast<int>((a * b) % kMod);
}

// Brute-force validator for n <= 4.
int64 brute_C(int n, int k_wanted) {
    if (n > 4) return -1;
    const int qsz = (1 << n) - 1;
    const int total_families = (1 << qsz) - 1;

    std::vector<int> subsets;
    subsets.reserve(qsz);
    for (int mask = 1; mask < (1 << n); ++mask) subsets.push_back(mask);

    int64 count = 0;
    for (int fam_mask = 1; fam_mask <= total_families; ++fam_mask) {
        std::vector<int> verts;
        for (int i = 0; i < qsz; ++i) {
            if (fam_mask & (1 << i)) verts.push_back(subsets[i]);
        }

        const int vcount = static_cast<int>(verts.size());
        std::vector<char> vis(vcount, 0);
        int comps = 0;
        for (int i = 0; i < vcount; ++i) {
            if (vis[i]) continue;
            ++comps;
            std::queue<int> q;
            q.push(i);
            vis[i] = 1;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (int v = 0; v < vcount; ++v) {
                    if (vis[v]) continue;
                    if ((verts[u] & verts[v]) != 0) {
                        vis[v] = 1;
                        q.push(v);
                    }
                }
            }
        }
        if (comps == k_wanted) ++count;
    }
    return count;
}

struct PolyMul {
    int N;
    int threads;

    PolyMul(int n, int t) : N(n), threads(std::max(1, t)) {}

    std::vector<int> mul(const std::vector<int>& a, const std::vector<int>& b) const {
        std::vector<int> c(N + 1, 0);
        if (threads <= 1 || N < 800) {
            for (int i = 0; i <= N; ++i) {
                __int128 sum = 0;
                for (int j = 0; j <= i; ++j) sum += (__int128)a[j] * b[i - j];
                c[i] = static_cast<int>(sum % kMod);
            }
            return c;
        }

        const long long total_ops = 1LL * (N + 1) * (N + 2) / 2;
        std::vector<int> boundary(threads + 1);
        boundary[0] = -1;
        boundary[threads] = N;
        int cur = -1;
        auto prefix_ops = [](long long i) -> long long {
            if (i < 0) return 0;
            return (i + 1) * (i + 2) / 2;
        };
        for (int t = 1; t < threads; ++t) {
            long long target = total_ops * t / threads;
            while (cur < N && prefix_ops(cur) < target) ++cur;
            boundary[t] = cur;
        }

        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (int t = 0; t < threads; ++t) {
            int L = boundary[t] + 1;
            int R = boundary[t + 1] + 1;
            pool.emplace_back([&, L, R]() {
                for (int i = L; i < R; ++i) {
                    __int128 sum = 0;
                    for (int j = 0; j <= i; ++j) sum += (__int128)a[j] * b[i - j];
                    c[i] = static_cast<int>(sum % kMod);
                }
            });
        }
        for (auto& th : pool) th.join();
        return c;
    }
};

std::vector<int> poly_pow(std::vector<int> base, int exp, const PolyMul& pm) {
    const int N = pm.N;
    std::vector<int> res(N + 1, 0);
    res[0] = 1;

    auto is_identity = [&](const std::vector<int>& p) {
        if (p[0] != 1) return false;
        for (int i = 1; i <= N; ++i) {
            if (p[i] != 0) return false;
        }
        return true;
    };

    while (exp > 0) {
        if (exp & 1) {
            if (is_identity(res)) res = base;
            else res = pm.mul(res, base);
        }
        exp >>= 1;
        if (exp == 0) break;
        base = pm.mul(base, base);
    }
    return res;
}

int solve(int N, int K, int threads, bool validate) {
    int max_n = std::max(N, K);
    if (validate) max_n = std::max(max_n, 100);

    std::vector<int> fact(max_n + 1), invfact(max_n + 1), invInt(max_n + 1);
    fact[0] = 1;
    for (int i = 1; i <= max_n; ++i) fact[i] = mulmod(fact[i - 1], i);
    invfact[max_n] = mod_inv(fact[max_n]);
    for (int i = max_n; i >= 1; --i) invfact[i - 1] = mulmod(invfact[i], i);
    invInt[1] = 1;
    for (int i = 2; i <= max_n; ++i) {
        invInt[i] = kMod - static_cast<int>(1LL * (kMod / i) * invInt[kMod % i] % kMod);
    }

    std::vector<int> pow2_modM1(max_n + 1, 1);
    for (int r = 1; r <= max_n; ++r) {
        pow2_modM1[r] = static_cast<int>(2LL * pow2_modM1[r - 1] % kModMinus1);
    }

    std::vector<int> T(max_n + 1, 0);
    for (int r = 0; r <= max_n; ++r) {
        int e = pow2_modM1[r] - 1;
        if (e < 0) e += kModMinus1;
        int64 v = mod_pow(2, e, kMod);
        int val = static_cast<int>(v) - 1;
        if (val < 0) val += kMod;
        T[r] = val;
    }

    if (validate && T[0] != 0) throw std::runtime_error("T(0) != 0");

    PolyMul pm(max_n, threads);

    std::vector<int> invfactSeries(max_n + 1), Bsign(max_n + 1);
    for (int i = 0; i <= max_n; ++i) invfactSeries[i] = invfact[i];
    for (int r = 0; r <= max_n; ++r) {
        int64 b = static_cast<int64>(T[r]) * invfact[r] % kMod;
        int bb = static_cast<int>(b);
        if (r & 1) bb = (bb == 0 ? 0 : kMod - bb);
        Bsign[r] = bb;
    }

    std::vector<int> convA = pm.mul(Bsign, invfactSeries);
    std::vector<int> Acoef(max_n + 1, 0);
    Acoef[0] = 1;  // empty structure for the exponential formula.
    for (int n = 1; n <= max_n; ++n) {
        int v = convA[n];
        if (n & 1) v = (v == 0 ? 0 : kMod - v);
        Acoef[n] = v;
    }

    std::vector<int> Gcoef(max_n + 1, 0), iG(max_n + 1, 0);
    for (int m = 1; m <= max_n; ++m) {
        int64 sum = 0;
        for (int i = 1; i < m; ++i) {
            sum += static_cast<int64>(iG[i]) * Acoef[m - i] % kMod;
            if (sum >= kMod) sum -= kMod;
        }
        int64 corr = sum * invInt[m] % kMod;
        int gm = static_cast<int>(Acoef[m] - corr);
        if (gm < 0) gm += kMod;
        Gcoef[m] = gm;
        iG[m] = static_cast<int>(1LL * m * Gcoef[m] % kMod);
    }

    if (validate && Gcoef[1] != 1) throw std::runtime_error("Gcoef[1] != 1");

    std::vector<int> Gpow = poly_pow(Gcoef, K, pm);
    std::vector<int> S = pm.mul(invfactSeries, Gpow);
    int answer = static_cast<int>(1LL * fact[N] * S[N] % kMod * invfact[K] % kMod);

    if (validate) {
        if (brute_C(2, 1) != 6) throw std::runtime_error("brute C(2,1) mismatch");
        if (brute_C(3, 1) != 111) throw std::runtime_error("brute C(3,1) mismatch");
        if (brute_C(4, 2) != 486) throw std::runtime_error("brute C(4,2) mismatch");
        if (N >= 100 && K == 10) {
            int chk = static_cast<int>(1LL * fact[100] * S[100] % kMod * invfact[10] % kMod);
            if (chk != 728209718) throw std::runtime_error("C(100,10) mismatch");
        }
        std::cerr << "Validation checkpoints passed.\n";
    }

    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int N = 10'000;
    int K = 10;
    unsigned hw = std::thread::hardware_concurrency();
    int threads = hw ? static_cast<int>(hw) : 1;
    threads = std::max(1, std::min(threads, 8));
    bool validate = true;

    // Optional CLI: ./a.out [N] [K] [threads] [validate(0/1)]
    if (argc >= 2) N = std::stoi(argv[1]);
    if (argc >= 3) K = std::stoi(argv[2]);
    if (argc >= 4) threads = std::max(1, std::stoi(argv[3]));
    if (argc >= 5) validate = (std::stoi(argv[4]) != 0);

    try {
        int answer = solve(N, K, threads, validate);
        std::cout << answer << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Validation/runtime error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
