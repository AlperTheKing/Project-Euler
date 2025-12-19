#include <bits/stdc++.h>
using namespace std;

static const long long MOD = 1000000007LL;
static const long long A = 1000000007LL;

struct PairData {
    int g = 1;
    long long mul = 0;
    vector<int> counts;
};

bool validate_permutation(const vector<int> &p, int n, const string &name) {
    vector<char> seen(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        int v = p[i];
        if (v < 1 || v > n) {
            cerr << "Validation failed: " << name << " has out-of-range value.\n";
            return false;
        }
        if (seen[v]) {
            cerr << "Validation failed: " << name << " has duplicate value.\n";
            return false;
        }
        seen[v] = 1;
    }
    return true;
}

vector<int> compute_counts(const vector<int> &Acyc, const vector<int> &Bcyc, int g) {
    const int s = static_cast<int>(Acyc.size());
    const int t = static_cast<int>(Bcyc.size());
    vector<vector<int>> A_r(g), B_r(g);
    for (int i = 0; i < s; ++i) {
        A_r[i % g].push_back(Acyc[i]);
    }
    for (int i = 0; i < t; ++i) {
        B_r[i % g].push_back(Bcyc[i]);
    }
    for (int r = 0; r < g; ++r) {
        sort(A_r[r].begin(), A_r[r].end());
        sort(B_r[r].begin(), B_r[r].end());
    }
    vector<vector<int>> f(g, vector<int>(g, 0));
    for (int r = 0; r < g; ++r) {
        for (int rp = 0; rp < g; ++rp) {
            int cnt = 0;
            const auto &a = A_r[r];
            const auto &b = B_r[rp];
            int j = 0;
            for (int x : a) {
                while (j < static_cast<int>(b.size()) && b[j] < x) {
                    ++j;
                }
                cnt += j;
            }
            f[r][rp] = cnt;
        }
    }
    // Count for any shift depends only on the offset modulo g.
    vector<int> counts(g, 0);
    for (int rd = 0; rd < g; ++rd) {
        int total = 0;
        for (int r = 0; r < g; ++r) {
            total += f[r][(r + rd) % g];
        }
        counts[rd] = total;
    }
    return counts;
}

int main() {
    const int m = 100;
    const int n = m * (m + 1) / 2;

    vector<int> tau(n + 1, 0), tau_inv(n + 1, 0), sigma(n + 1, 0), pi(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        tau[i] = static_cast<int>((A * i) % n) + 1;
    }
    for (int i = 1; i <= n; ++i) {
        tau_inv[tau[i]] = i;
    }
    for (int i = 1; i <= n; ++i) {
        if (tau_inv[tau[i]] != i) {
            cerr << "Validation failed: tau inverse mismatch.\n";
            return 0;
        }
    }

    for (int i = 1; i < n; ++i) {
        sigma[i] = i + 1;
    }
    int prev = 0;
    for (int k = 1; k <= m; ++k) {
        int t = k * (k + 1) / 2;
        sigma[t] = prev + 1;
        prev = t;
    }

    for (int i = 1; i <= n; ++i) {
        pi[i] = tau_inv[sigma[tau[i]]];
    }
    if (!validate_permutation(pi, n, "pi")) {
        return 0;
    }

    vector<vector<int>> cycles;
    vector<int> cycle_id(n + 1, -1), pos_in(n + 1, -1);
    vector<char> visited(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        if (visited[i]) {
            continue;
        }
        vector<int> cyc;
        int cur = i;
        while (!visited[cur]) {
            visited[cur] = 1;
            pos_in[cur] = static_cast<int>(cyc.size());
            cycle_id[cur] = static_cast<int>(cycles.size());
            cyc.push_back(cur);
            cur = pi[cur];
        }
        cycles.push_back(move(cyc));
    }

    int total_len = 0;
    vector<int> len_count(m + 1, 0);
    for (const auto &cyc : cycles) {
        total_len += static_cast<int>(cyc.size());
        if (static_cast<int>(cyc.size()) >= 1 && static_cast<int>(cyc.size()) <= m) {
            len_count[cyc.size()]++;
        }
    }
    if (total_len != n) {
        cerr << "Validation failed: cycle coverage mismatch.\n";
        return 0;
    }
    for (int len = 1; len <= m; ++len) {
        if (len_count[len] != 1) {
            cerr << "Validation failed: unexpected cycle length distribution.\n";
            return 0;
        }
    }

    vector<long long> fact(n + 1, 1), weight(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        fact[i] = (fact[i - 1] * i) % MOD;
    }
    for (int i = 1; i <= n; ++i) {
        weight[i] = fact[n - i];
    }

    long long m_fact = 1;
    for (int i = 1; i <= m; ++i) {
        m_fact = (m_fact * i) % MOD;
    }

    const int num_cycles = static_cast<int>(cycles.size());
    vector<int> cycle_len(num_cycles, 0);
    for (int i = 0; i < num_cycles; ++i) {
        cycle_len[i] = static_cast<int>(cycles[i].size());
    }

    int maxL = 1;
    for (int i = 0; i < num_cycles; ++i) {
        for (int j = 0; j < num_cycles; ++j) {
            int g = std::gcd(cycle_len[i], cycle_len[j]);
            int L = cycle_len[i] / g * cycle_len[j];
            maxL = max(maxL, L);
        }
    }
    vector<long long> inv(maxL + 1, 1);
    for (int i = 2; i <= maxL; ++i) {
        inv[i] = MOD - (MOD / i) * inv[MOD % i] % MOD;
    }

    vector<vector<PairData>> pair_data(num_cycles, vector<PairData>(num_cycles));
    for (int ci = 0; ci < num_cycles; ++ci) {
        for (int cj = 0; cj < num_cycles; ++cj) {
            int s = cycle_len[ci];
            int t = cycle_len[cj];
            int g = std::gcd(s, t);
            int L = s / g * t;
            PairData pd;
            pd.g = g;
            pd.mul = (m_fact * inv[L]) % MOD;
            pd.counts = compute_counts(cycles[ci], cycles[cj], g);
            pair_data[ci][cj] = move(pd);
        }
    }

    long long total = m_fact % MOD;
    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 1;
    }
    if (threads > 8) {
        threads = 8;
    }
    vector<long long> partial(threads, 0);
    vector<thread> workers;

    int block = n / static_cast<int>(threads);
    int start = 1;
    for (unsigned int t = 0; t < threads; ++t) {
        int end = (t == threads - 1) ? n : (start + block - 1);
        workers.emplace_back([&, start, end, t]() {
            long long local = 0;
            for (int i = start; i <= end; ++i) {
                long long wi = weight[i];
                int ci = cycle_id[i];
                int pos_i = pos_in[i];
                for (int j = i + 1; j <= n; ++j) {
                    int cj = cycle_id[j];
                    int pos_j = pos_in[j];
                    const PairData &pd = pair_data[ci][cj];
                    int g = pd.g;
                    int r = pos_j - pos_i;
                    r %= g;
                    if (r < 0) {
                        r += g;
                    }
                    int cnt = pd.counts[r];
                    long long add = (wi * (static_cast<long long>(cnt) * pd.mul % MOD)) % MOD;
                    local += add;
                    if (local >= MOD) {
                        local -= MOD;
                    }
                }
            }
            partial[t] = local;
        });
        start = end + 1;
    }
    for (auto &th : workers) {
        th.join();
    }
    for (long long v : partial) {
        total += v;
        if (total >= MOD) {
            total -= MOD;
        }
    }

    cout << (total % MOD) << "\n";
    return 0;
}
