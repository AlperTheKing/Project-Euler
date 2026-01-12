#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <thread>
#include <vector>

struct Transition {
    int next;
    int add;
    long double weight;
};

class Solver {
public:
    Solver(int days, int cluster_size, int within_days)
        : N(days), k(cluster_size), D(within_days), W(within_days + 1) {
        L = W - 1;
        n_max = (k - 1) * N / W;
        build_inv_factorials();
        build_states();
        build_transitions();
    }

    long double solve() const {
        const int S = static_cast<int>(states.size());
        unsigned threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 4;
        if (threads > static_cast<unsigned>(S)) threads = S;

        std::vector<std::vector<long double>> partial(threads, std::vector<long double>(n_max + 1, 0.0L));
        std::vector<std::thread> pool;
        pool.reserve(threads);

        auto worker = [&](unsigned tid, int start, int end) {
            std::vector<long double> dp(S * (n_max + 1), 0.0L);
            std::vector<long double> next(S * (n_max + 1), 0.0L);

            for (int s = start; s < end; ++s) {
                std::fill(dp.begin(), dp.end(), 0.0L);
                dp[s * (n_max + 1)] = 1.0L;

                for (int day = 0; day < N; ++day) {
                    std::fill(next.begin(), next.end(), 0.0L);
                    for (int u = 0; u < S; ++u) {
                        const long double* row = &dp[u * (n_max + 1)];
                        for (const auto& tr : transitions[u]) {
                            long double* dest = &next[tr.next * (n_max + 1)];
                            int limit = n_max - tr.add;
                            for (int n = 0; n <= limit; ++n) {
                                long double val = row[n];
                                if (val != 0.0L) {
                                    dest[n + tr.add] += val * tr.weight;
                                }
                            }
                        }
                    }
                    dp.swap(next);
                }

                const long double* row = &dp[s * (n_max + 1)];
                for (int n = 0; n <= n_max; ++n) {
                    partial[tid][n] += row[n];
                }
            }
        };

        int chunk = (S + static_cast<int>(threads) - 1) / static_cast<int>(threads);
        for (unsigned t = 0; t < threads; ++t) {
            int start = static_cast<int>(t) * chunk;
            int end = std::min(start + chunk, S);
            if (start >= end) break;
            pool.emplace_back(worker, t, start, end);
        }
        for (auto& th : pool) th.join();

        std::vector<long double> totals(n_max + 1, 0.0L);
        for (unsigned t = 0; t < threads; ++t) {
            for (int n = 0; n <= n_max; ++n) {
                totals[n] += partial[t][n];
            }
        }

        long double expected = 0.0L;
        long double fact = 1.0L;
        long double invN = 1.0L / static_cast<long double>(N);
        long double invNpow = 1.0L;
        for (int n = 0; n <= n_max; ++n) {
            if (n > 0) {
                fact *= static_cast<long double>(n);
                invNpow *= invN;
            }
            long double prob = totals[n] * fact * invNpow;
            expected += prob;
        }
        return expected;
    }

private:
    int N = 0;
    int k = 0;
    int D = 0;
    int W = 0;
    int L = 0;
    int n_max = 0;
    std::vector<std::vector<int>> states;
    std::map<std::vector<int>, int> state_id;
    std::vector<std::vector<Transition>> transitions;
    std::vector<long double> inv_fact;

    void build_inv_factorials() {
        inv_fact.assign(k, 1.0L);
        long double fact = 1.0L;
        for (int i = 1; i < k; ++i) {
            fact *= static_cast<long double>(i);
            inv_fact[i] = 1.0L / fact;
        }
    }

    void build_states() {
        states.clear();
        state_id.clear();
        if (L == 0) {
            states.emplace_back();
            state_id[states.back()] = 0;
            return;
        }

        std::vector<int> cur(L, 0);
        generate_states(0, k - 1, cur);
    }

    void generate_states(int pos, int remaining, std::vector<int>& cur) {
        if (pos == L) {
            int id = static_cast<int>(states.size());
            states.push_back(cur);
            state_id[cur] = id;
            return;
        }
        for (int v = 0; v <= remaining; ++v) {
            cur[pos] = v;
            generate_states(pos + 1, remaining - v, cur);
        }
    }

    void build_transitions() {
        transitions.assign(states.size(), {});
        for (size_t i = 0; i < states.size(); ++i) {
            int sum = 0;
            for (int x : states[i]) sum += x;
            for (int c = 0; c + sum <= k - 1; ++c) {
                std::vector<int> next;
                if (L > 0) {
                    next.assign(states[i].begin() + 1, states[i].end());
                    next.push_back(c);
                }
                auto it = state_id.find(next);
                transitions[i].push_back({it->second, c, inv_fact[c]});
            }
        }
    }
};

static void run_check(const char* name, int days, int cluster, int within, long double expected) {
    Solver solver(days, cluster, within);
    long double result = solver.solve();
    long double diff = std::fabsl(result - expected);
    std::cout << "[" << name << "] "
              << "N=" << days << " k=" << cluster << " within=" << within
              << " -> " << std::fixed << std::setprecision(8) << result
              << " (expected " << std::fixed << std::setprecision(8) << expected << ")";
    if (diff < 1e-7L) {
        std::cout << " [PASS]\n";
    } else {
        std::cout << " [FAIL]\n";
    }
}

int main() {
    run_check("WimWi", 10, 3, 1, 5.78688636L);
    run_check("Joka", 100, 3, 7, 8.48967364L);

    std::cout << "------------------------------------------------\n";
    Solver earth(365, 4, 7);
    long double ans = earth.solve();
    std::cout << "Earth answer: " << std::fixed << std::setprecision(8) << ans << "\n";
    return 0;
}
