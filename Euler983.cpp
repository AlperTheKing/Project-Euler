#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

using namespace std;
using u64 = unsigned long long;
using i128 = __int128_t;

static inline u64 mul_cap(u64 a, u64 b, u64 cap) {
    if (a == 0 || b == 0) return 0;
    if (a > cap / b) return cap;
    return a * b;
}

static vector<int> primes_1mod4(int limit) {
    vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; 1LL * i * i <= limit; i++) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= limit; j += i) is_prime[j] = false;
    }
    vector<int> ps;
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i] && i % 4 == 1) ps.push_back(i);
    }
    return ps;
}

static inline u64 get_best_atomic(const atomic<u64>& best) {
    return best.load(memory_order_relaxed);
}

static void update_best_atomic(atomic<u64>& best, u64 cand) {
    u64 cur = best.load(memory_order_relaxed);
    while (cand < cur && !best.compare_exchange_weak(cur, cand, memory_order_relaxed)) {
    }
}

// Minimize product p_i^e_i with non-increasing exponents and Π(e_i+1) >= target.
static void dfs_minimize(const vector<int>& ps,
                         int idx,
                         int max_exp,
                         i128 val,
                         u64 div_prod,
                         u64 target_div_prod,
                         atomic<u64>& best) {
    if (div_prod >= target_div_prod) {
        update_best_atomic(best, static_cast<u64>(val));
        return;
    }

    if (val >= static_cast<i128>(get_best_atomic(best))) return;
    if (idx >= static_cast<int>(ps.size())) return;

    int rem = static_cast<int>(ps.size()) - idx;
    u64 max_possible = div_prod;
    for (int i = 0; i < rem && max_possible < target_div_prod; i++) {
        max_possible = mul_cap(max_possible, static_cast<u64>(max_exp + 1), target_div_prod);
    }
    if (max_possible < target_div_prod) return;

    i128 cur = val;
    for (int e = 1; e <= max_exp; e++) {
        cur *= static_cast<i128>(ps[idx]);
        if (cur >= static_cast<i128>(get_best_atomic(best))) break;
        dfs_minimize(ps, idx + 1, e, cur, div_prod * static_cast<u64>(e + 1), target_div_prod, best);
    }
}

int main() {
    const int n = 500;
    // Need r2(m) = 4*Π(e+1) >= n.
    const u64 target_div_prod = (n + 3) / 4;  // ceil(n/4)

    auto ps = primes_1mod4(3000);
    assert(!ps.empty() && ps[0] == 5);

    // Initial bound: exponents 1 on first k primes.
    u64 init_div_prod = 1;
    i128 init_val = 1;
    int k = 0;
    while (init_div_prod < target_div_prod) {
        init_val *= static_cast<i128>(ps[k]);
        init_div_prod *= 2;
        k++;
    }
    assert(init_div_prod >= target_div_prod);

    atomic<u64> best;
    best.store(static_cast<u64>(init_val), memory_order_relaxed);

    int max_e0 = 0;
    {
        i128 v = 1;
        while (true) {
            i128 nv = v * static_cast<i128>(5);
            if (nv >= static_cast<i128>(get_best_atomic(best))) break;
            v = nv;
            max_e0++;
        }
    }
    max_e0 = max(1, max_e0);

    unsigned hw = thread::hardware_concurrency();
    if (hw == 0) hw = 4;
    mutex qmtx;
    int next_e = 1;

    auto worker = [&]() {
        while (true) {
            int e0;
            {
                lock_guard<mutex> lock(qmtx);
                if (next_e > max_e0) return;
                e0 = next_e++;
            }
            i128 val = 1;
            for (int i = 0; i < e0; i++) val *= static_cast<i128>(5);
            if (val >= static_cast<i128>(get_best_atomic(best))) continue;
            dfs_minimize(ps, 1, e0, val, static_cast<u64>(e0 + 1), target_div_prod, best);
        }
    };

    vector<thread> threads;
    threads.reserve(hw);
    for (unsigned i = 0; i < hw; i++) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    u64 ans = best.load(memory_order_relaxed);

    // Validation checkpoints.
    u64 tmp = ans;
    u64 div_check = 1;
    for (int p : ps) {
        u64 pp = static_cast<u64>(p);
        if (pp * pp > tmp) break;
        if (tmp % pp == 0) {
            int e = 0;
            while (tmp % pp == 0) {
                tmp /= pp;
                e++;
            }
            div_check *= static_cast<u64>(e + 1);
        }
    }
    if (tmp > 1) div_check *= 2;
    assert(div_check >= target_div_prod);
    assert(4 * div_check >= static_cast<u64>(n));

    cout << ans << "\n";
    return 0;
}
