#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

static constexpr double EPS = 1e-10;

template <typename Func>
static void parallel_for(int total, int threads, Func func) {
    if (threads <= 1 || total < 128) {
        for (int i = 0; i < total; ++i) func(i);
        return;
    }
    threads = min(threads, total);
    int chunk = (total + threads - 1) / threads;
    vector<thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t) {
        int start = t * chunk;
        int end = min(total, start + chunk);
        if (start >= end) break;
        pool.emplace_back([=, &func]() {
            for (int i = start; i < end; ++i) func(i);
        });
    }
    for (auto &th : pool) th.join();
}

struct Simplex {
    int m;
    int n;
    vector<vector<double>> T;
    vector<int> basis;

    Simplex(int rows, int cols)
        : m(rows), n(cols), T(rows + 1, vector<double>(cols + 1, 0.0)), basis(rows, -1) {}

    void pivot(int r, int c) {
        double inv = 1.0 / T[r][c];
        for (int j = 0; j <= n; ++j) T[r][j] *= inv;
        for (int i = 0; i <= m; ++i) {
            if (i == r) continue;
            double factor = T[i][c];
            if (fabs(factor) <= EPS) continue;
            for (int j = 0; j <= n; ++j) T[i][j] -= factor * T[r][j];
        }
        basis[r] = c;
    }

    void set_objective(const vector<double> &c) {
        fill(T[m].begin(), T[m].end(), 0.0);
        for (int j = 0; j < n; ++j) T[m][j] = -c[j];
        T[m][n] = 0.0;
        for (int i = 0; i < m; ++i) {
            int var = basis[i];
            if (var < 0) continue;
            double coeff = c[var];
            if (fabs(coeff) <= EPS) continue;
            for (int j = 0; j <= n; ++j) T[m][j] += coeff * T[i][j];
        }
    }

    bool solve() {
        while (true) {
            int enter = -1;
            for (int j = 0; j < n; ++j) {
                if (T[m][j] < -EPS) {
                    enter = j;
                    break;
                }
            }
            if (enter == -1) return true;

            double min_ratio = 0.0;
            int leave = -1;
            for (int i = 0; i < m; ++i) {
                double a = T[i][enter];
                if (a > EPS) {
                    double ratio = T[i][n] / a;
                    if (leave == -1 || ratio < min_ratio - 1e-12 ||
                        (fabs(ratio - min_ratio) <= 1e-12 && basis[i] < basis[leave])) {
                        min_ratio = ratio;
                        leave = i;
                    }
                }
            }
            if (leave == -1) return false;
            pivot(leave, enter);
        }
    }

    double objective_value() const { return T[m][n]; }
};

struct GameData {
    int dice = 0;
    int actions = 0;
    int states = 0;
    int obs = 0;
    double prob = 0.0;
    vector<int> visible_max;
    vector<int> action_obs;
    vector<int> action_hidden;
};

static GameData build_game(int dice) {
    GameData g;
    g.dice = dice;
    g.actions = dice;
    if (dice == 2) {
        g.states = 36;
        g.obs = 6;
        g.prob = 1.0 / 36.0;
        g.visible_max.resize(g.obs);
        for (int v = 0; v < g.obs; ++v) g.visible_max[v] = v + 1;

        int n_action = g.states * g.actions;
        g.action_obs.resize(n_action);
        g.action_hidden.resize(n_action);

        for (int a = 1; a <= 6; ++a) {
            for (int b = 1; b <= 6; ++b) {
                int s = (a - 1) * 6 + (b - 1);
                int base = s * 2;
                g.action_obs[base] = b - 1;
                g.action_hidden[base] = a;
                g.action_obs[base + 1] = a - 1;
                g.action_hidden[base + 1] = b;
            }
        }
        return g;
    }
    if (dice == 3) {
        g.states = 216;
        g.obs = 21;
        g.prob = 1.0 / 216.0;

        int obs_index[7][7];
        for (int i = 0; i <= 6; ++i) for (int j = 0; j <= 6; ++j) obs_index[i][j] = -1;

        g.visible_max.resize(g.obs);
        int idx = 0;
        for (int u = 1; u <= 6; ++u) {
            for (int v = u; v <= 6; ++v) {
                obs_index[u][v] = idx;
                obs_index[v][u] = idx;
                g.visible_max[idx] = v;
                idx++;
            }
        }

        int n_action = g.states * g.actions;
        g.action_obs.resize(n_action);
        g.action_hidden.resize(n_action);

        for (int a = 1; a <= 6; ++a) {
            for (int b = 1; b <= 6; ++b) {
                for (int c = 1; c <= 6; ++c) {
                    int s = ((a - 1) * 6 + (b - 1)) * 6 + (c - 1);
                    int base = s * 3;
                    g.action_obs[base] = obs_index[b][c];
                    g.action_hidden[base] = a;
                    g.action_obs[base + 1] = obs_index[a][c];
                    g.action_hidden[base + 1] = b;
                    g.action_obs[base + 2] = obs_index[a][b];
                    g.action_hidden[base + 2] = c;
                }
            }
        }
        return g;
    }
    cerr << "Unsupported dice count: " << dice << "\n";
    return g;
}

struct LPData {
    Simplex lp;
    int n_action;
    int n_T;
    int n_slack;
    int n_art;
    int art_start;
    int slack_start;

    LPData(int rows, int cols, int actions, int tvars, int slack, int art, int art_s, int slack_s)
        : lp(rows, cols), n_action(actions), n_T(tvars), n_slack(slack), n_art(art),
          art_start(art_s), slack_start(slack_s) {}
};

static LPData build_lp(const GameData &g, unsigned threads) {
    int n_action = g.states * g.actions;
    int n_T = g.obs;
    int n_orig = n_action + n_T;
    int n_slack = g.obs * 2;
    int n_art = g.states;
    int n_total = n_orig + n_slack + n_art;
    int rows = g.states + n_slack;

    int slack_start = n_orig;
    int art_start = n_orig + n_slack;

    LPData data(rows, n_total, n_action, n_T, n_slack, n_art, art_start, slack_start);
    Simplex &lp = data.lp;

    for (int s = 0; s < g.states; ++s) {
        int row = s;
        int base = s * g.actions;
        for (int a = 0; a < g.actions; ++a) lp.T[row][base + a] = 1.0;
        lp.T[row][art_start + s] = 1.0;
        lp.T[row][n_total] = 1.0;
        lp.basis[row] = art_start + s;
    }

    for (int o = 0; o < g.obs; ++o) {
        for (int t = 0; t < 2; ++t) {
            int row = g.states + o * 2 + t;
            int tvar = n_action + o;
            lp.T[row][tvar] = -1.0;
            lp.T[row][slack_start + o * 2 + t] = 1.0;
            lp.T[row][n_total] = 0.0;
            lp.basis[row] = slack_start + o * 2 + t;
        }
    }

    auto worker = [&](int idx) {
        int obs = g.action_obs[idx];
        int hidden = g.action_hidden[idx];
        int row_hidden = g.states + obs * 2;
        int row_visible = row_hidden + 1;
        lp.T[row_hidden][idx] += g.prob * hidden;
        lp.T[row_visible][idx] += g.prob * g.visible_max[obs];
    };

    parallel_for(n_action, static_cast<int>(threads), worker);

    return data;
}

static void remove_artificial(Simplex &lp, int art_start, int art_count) {
    if (art_count == 0) return;
    vector<int> rows_to_remove;
    vector<char> is_art(lp.n, 0);
    for (int i = 0; i < art_count; ++i) is_art[art_start + i] = 1;

    for (int i = 0; i < lp.m; ++i) {
        int var = lp.basis[i];
        if (var < 0 || !is_art[var]) continue;
        int pivot_col = -1;
        for (int j = 0; j < lp.n; ++j) {
            if (is_art[j]) continue;
            if (fabs(lp.T[i][j]) > EPS) {
                pivot_col = j;
                break;
            }
        }
        if (pivot_col != -1) {
            lp.pivot(i, pivot_col);
        } else {
            rows_to_remove.push_back(i);
        }
    }

    if (!rows_to_remove.empty()) {
        vector<char> drop(lp.m, 0);
        for (int r : rows_to_remove) drop[r] = 1;
        int new_m = lp.m - static_cast<int>(rows_to_remove.size());
        vector<vector<double>> newT(new_m + 1, vector<double>(lp.n + 1, 0.0));
        vector<int> new_basis(new_m, -1);
        int r2 = 0;
        for (int i = 0; i < lp.m; ++i) {
            if (drop[i]) continue;
            newT[r2] = lp.T[i];
            new_basis[r2] = lp.basis[i];
            r2++;
        }
        newT[new_m] = lp.T[lp.m];
        lp.T.swap(newT);
        lp.basis.swap(new_basis);
        lp.m = new_m;
    }

    int new_n = lp.n - art_count;
    vector<vector<double>> newT(lp.m + 1, vector<double>(new_n + 1, 0.0));
    for (int i = 0; i <= lp.m; ++i) {
        int col_new = 0;
        for (int j = 0; j < lp.n; ++j) {
            if (is_art[j]) continue;
            newT[i][col_new++] = lp.T[i][j];
        }
        newT[i][new_n] = lp.T[i][lp.n];
    }

    lp.n = new_n;
    lp.T.swap(newT);
}

static double solve_game(int dice, unsigned threads) {
    GameData g = build_game(dice);
    if (g.states == 0) return nan("");

    LPData data = build_lp(g, threads);
    Simplex &lp = data.lp;

    vector<double> c1(lp.n, 0.0);
    for (int i = 0; i < data.n_art; ++i) c1[data.art_start + i] = -1.0;
    lp.set_objective(c1);
    if (!lp.solve()) {
        cerr << "Phase I unbounded\n";
        return nan("");
    }

    double phase1 = lp.objective_value();
    if (phase1 < -1e-8) {
        cerr << "Phase I infeasible: " << phase1 << "\n";
        return nan("");
    }

    remove_artificial(lp, data.art_start, data.n_art);

    vector<double> c2(lp.n, 0.0);
    for (int i = 0; i < data.n_T; ++i) c2[data.n_action + i] = -1.0;
    lp.set_objective(c2);
    if (!lp.solve()) {
        cerr << "Phase II unbounded\n";
        return nan("");
    }

    return -lp.objective_value();
}

static bool run_validation(unsigned threads) {
    double expected = 145.0 / 36.0;
    double got = solve_game(2, threads);
    if (!isfinite(got) || fabs(got - expected) > 1e-7) {
        cerr << "Validation failed: expected " << fixed << setprecision(9) << expected
             << ", got " << got << "\n";
        return false;
    }
    return true;
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    unsigned threads = thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = min(threads, 8u);
    bool validate = true;

    if (argc >= 2) threads = max(1u, static_cast<unsigned>(stoul(argv[1])));
    if (argc >= 3) validate = (stoul(argv[2]) != 0);

    if (validate && !run_validation(min(threads, 4u))) {
        return 1;
    }

    double answer = solve_game(3, threads);
    if (!isfinite(answer)) return 1;

    cout << fixed << setprecision(6) << answer << "\n";
    return 0;
}
