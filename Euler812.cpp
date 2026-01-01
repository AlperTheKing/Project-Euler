#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

static constexpr int MOD = 998244353;

static inline int addmod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static inline int submod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

static vector<int> compute_phi(int M) {
    vector<int> phi(M + 1);
    for (int i = 0; i <= M; ++i) phi[i] = i;
    for (int i = 2; i <= M; ++i) {
        if (phi[i] == i) {
            for (int j = i; j <= M; j += i) {
                phi[j] -= phi[j] / i;
            }
        }
    }
    return phi;
}

// Chain for odd base b=1 has the special constraint a_4 <= 2*a_2.
// This function builds its full generating polynomial up to degree n.
static vector<int> build_chain1(int n, const vector<int>& deg) {
    vector<int> tail_weights;
    long long prefix = 0;
    for (long long m = 4; m < static_cast<long long>(deg.size()); m *= 2) {
        int d = deg[static_cast<size_t>(m)];
        if (d > n) break;
        prefix += d;
        long long t = prefix + 1; // y=x substitution contributes one extra degree per count.
        if (t > n) break;
        tail_weights.push_back(static_cast<int>(t));
    }

    vector<int> G(n + 1, 0);
    vector<int> H(n + 1, 0);
    G[0] = 1;
    H[0] = 1;

    for (int t : tail_weights) {
        for (int k = t; k <= n; ++k) {
            G[k] = addmod(G[k], G[k - t]);
        }
    }

    for (int t : tail_weights) {
        for (int k = t; k <= n; ++k) {
            H[k] = submod(H[k], H[k - t]);
        }
    }

    const int inv2 = (MOD + 1) / 2;
    vector<int> S(n + 1, 0);
    for (int k = 0; k <= n; ++k) {
        long long val = static_cast<long long>(G[k]) + H[k];
        if (k > 0) {
            val += static_cast<long long>(G[k - 1]) - H[k - 1];
        }
        val %= MOD;
        if (val < 0) val += MOD;
        S[k] = static_cast<int>(val * inv2 % MOD);
    }

    for (int k = 1; k <= n; ++k) {
        S[k] = addmod(S[k], S[k - 1]);
    }
    for (int k = 2; k <= n; ++k) {
        S[k] = addmod(S[k], S[k - 2]);
    }

    return S;
}

static int compute_S(int n) {
    const int M = max(100, 40 * n); // Safe envelope for phi(m) <= 2n in this problem range.
    const vector<int> phi = compute_phi(M);
    vector<int> deg(M + 1, 0);
    if (M >= 1) deg[1] = 1;
    if (M >= 2) deg[2] = 1;
    for (int m = 3; m <= M; ++m) {
        deg[m] = phi[m] / 2;
    }

    vector<int> dp = build_chain1(n, deg);

    for (int b = 3; b <= M; b += 2) {
        if (deg[b] > n) continue;
        long long prefix = 0;
        long long m = b;
        while (m <= M) {
            int d = deg[static_cast<size_t>(m)];
            if (d > n) break;
            prefix += d;
            if (prefix > n) break;
            const int w = static_cast<int>(prefix);
            for (int k = w; k <= n; ++k) {
                dp[k] = addmod(dp[k], dp[k - w]);
            }
            m *= 2;
        }
    }

    return dp[n];
}

static bool run_validation() {
    bool ok = true;
    if (compute_S(2) != 6) {
        cerr << "Validation failed: S(2) != 6\n";
        ok = false;
    }
    if (compute_S(5) != 58) {
        cerr << "Validation failed: S(5) != 58\n";
        ok = false;
    }
    if (compute_S(20) != 122087) {
        cerr << "Validation failed: S(20) != 122087\n";
        ok = false;
    }
    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 10000;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--validate") {
            validate = true;
        } else {
            n = stoi(arg);
        }
    }

    if (validate && !run_validation()) {
        return 1;
    }

    cout << compute_S(n) << "\n";
    return 0;
}
