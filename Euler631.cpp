#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr uint64_t MOD = 1000000007ULL;
constexpr int MAX_PERM = 64;

struct State {
    int n = 0;
    int inv = 0;
    array<int, MAX_PERM> perm{};
};

template <typename F>
void for_each_child(State& st, int m, int max_len, F&& func) {
    const int n = st.n;
    if (n >= max_len) return;

    const int INF = 1'000'000'000;
    array<int, MAX_PERM> min_nonrecord{};
    array<int, MAX_PERM> max_suffix{};

    // min_nonrecord[p] tracks the smallest prefix element that has a smaller predecessor.
    int min_val = INF;
    int min_nr = INF;
    min_nonrecord[0] = INF;
    for (int i = 0; i < n; ++i) {
        int v = st.perm[i];
        if (v > min_val) {
            if (v < min_nr) min_nr = v;
        } else {
            min_val = v;
        }
        min_nonrecord[i + 1] = min_nr;
    }

    int max_val = -1;
    max_suffix[n] = -1;
    for (int i = n - 1; i >= 0; --i) {
        int v = st.perm[i];
        if (v > max_val) max_val = v;
        max_suffix[i] = max_val;
    }

    const int val = n + 1;
    for (int pos = 0; pos <= n; ++pos) {
        int mn = min_nonrecord[pos];
        // Inserting the new max before pos is safe iff suffix max <= min_nonrecord(prefix).
        if (mn != INF && max_suffix[pos] > mn) continue;
        int add = n - pos;
        if (st.inv + add > m) continue;

        for (int i = n; i > pos; --i) st.perm[i] = st.perm[i - 1];
        st.perm[pos] = val;
        st.n += 1;
        st.inv += add;
        func();
        st.n -= 1;
        st.inv -= add;
        for (int i = pos; i < n; ++i) st.perm[i] = st.perm[i + 1];
    }
}

void dfs(State& st, int m, int max_len, vector<uint64_t>& counts) {
    counts[st.n] += 1;
    if (st.n >= max_len) return;
    for_each_child(st, m, max_len, [&]() { dfs(st, m, max_len, counts); });
}

void collect_tasks(State& st, int m, int max_len, int cut,
                   vector<uint64_t>& base_counts, vector<State>& tasks) {
    if (st.n < cut) base_counts[st.n] += 1;
    if (st.n >= max_len) return;
    if (st.n == cut) {
        tasks.push_back(st);
        return;
    }
    for_each_child(st, m, max_len,
                   [&]() { collect_tasks(st, m, max_len, cut, base_counts, tasks); });
}

vector<uint64_t> compute_counts(int m, int max_len, unsigned threads) {
    vector<uint64_t> counts(static_cast<size_t>(max_len + 1), 0);
    if (max_len == 0) {
        counts[0] = 1;
        return counts;
    }

    if (threads <= 1 || max_len <= 5) {
        State st;
        dfs(st, m, max_len, counts);
        return counts;
    }

    int cut = min(6, max_len);
    vector<uint64_t> base_counts(static_cast<size_t>(max_len + 1), 0);
    vector<State> tasks;
    State st;
    collect_tasks(st, m, max_len, cut, base_counts, tasks);

    if (tasks.empty()) return base_counts;

    unsigned tcount = min<unsigned>(threads, tasks.size());
    vector<vector<uint64_t>> thread_counts(
        tcount, vector<uint64_t>(static_cast<size_t>(max_len + 1), 0));
    atomic<size_t> task_index{0};
    vector<thread> workers;
    workers.reserve(tcount);

    for (unsigned tid = 0; tid < tcount; ++tid) {
        workers.emplace_back([&, tid]() {
            while (true) {
                size_t idx = task_index.fetch_add(1, memory_order_relaxed);
                if (idx >= tasks.size()) break;
                State local = tasks[idx];
                dfs(local, m, max_len, thread_counts[tid]);
            }
        });
    }
    for (auto& th : workers) th.join();

    counts = base_counts;
    for (unsigned tid = 0; tid < tcount; ++tid) {
        for (int i = 0; i <= max_len; ++i) counts[i] += thread_counts[tid][i];
    }
    return counts;
}

uint64_t sum_counts(const vector<uint64_t>& counts, int up_to) {
    uint64_t sum = 0;
    for (int i = 0; i <= up_to; ++i) sum += counts[i];
    return sum;
}

uint64_t sum_counts_mod(const vector<uint64_t>& counts, int up_to) {
    uint64_t sum = 0;
    for (int i = 0; i <= up_to; ++i) sum = (sum + counts[i]) % MOD;
    return sum;
}

uint64_t compute_f(uint64_t n, int m, unsigned threads, bool require_stable) {
    const int stable_len = m + 2;
    const int max_len = stable_len + 1;
    vector<uint64_t> counts = compute_counts(m, max_len, threads);

    // Once counts stabilize, appending the new maximum is the only valid extension.
    if (require_stable && counts[stable_len] != counts[stable_len + 1]) {
        cerr << "Stabilization check failed at length " << stable_len << ".\n";
        exit(1);
    }

    if (n < static_cast<uint64_t>(stable_len)) {
        return sum_counts_mod(counts, static_cast<int>(n));
    }

    uint64_t prefix = sum_counts_mod(counts, stable_len - 1);
    uint64_t g = counts[stable_len] % MOD;
    uint64_t mult = (n - static_cast<uint64_t>(stable_len - 1)) % MOD;
    return (prefix + (g * mult) % MOD) % MOD;
}

void run_checks() {
    struct Check {
        uint64_t n;
        int m;
        uint64_t expected;
    };

    vector<Check> checks = {
        {2, 0, 3},
        {4, 5, 32},
        {10, 25, 294400},
    };

    for (const auto& chk : checks) {
        int max_len = static_cast<int>(chk.n);
        vector<uint64_t> counts = compute_counts(chk.m, max_len, 1);
        uint64_t got = sum_counts(counts, max_len);
        if (got != chk.expected) {
            cerr << "Validation failed for f(" << chk.n << "," << chk.m << "). "
                 << "Expected " << chk.expected << ", got " << got << ".\n";
            exit(1);
        }
    }
}

}  // namespace

int main() {
    run_checks();

    const uint64_t n = 1000000000000000000ULL;
    const int m = 40;
    unsigned threads = thread::hardware_concurrency();
    if (threads == 0) threads = 1;

    uint64_t answer = compute_f(n, m, threads, true);
    cout << answer << "\n";
    return 0;
}
