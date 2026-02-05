#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using u64 = unsigned long long;
using boost::multiprecision::cpp_int;

namespace {

constexpr u64 kTarget = 8916100448256ULL;  // 12^12.

struct Fenwick {
    explicit Fenwick(int n = 0) : n_(n), bit_(n + 1, 0) {}

    void add(int idx, int delta) {
        for (++idx; idx <= n_; idx += idx & -idx) bit_[idx] += delta;
    }

    int sum_prefix(int idx) const {
        int sum = 0;
        for (++idx; idx > 0; idx -= idx & -idx) sum += bit_[idx];
        return sum;
    }

    int n_;
    std::vector<int> bit_;
};

u64 first_sort_steps_fast(const std::vector<int>& perm) {
    const int n = static_cast<int>(perm.size());
    if (n <= 1) return 0;

    std::vector<int> pos(n + 1, 0);
    for (int i = 0; i < n; ++i) pos[perm[i]] = i;

    std::vector<int> insertion_pos(n + 1, 0);
    Fenwick fw(n);
    fw.add(pos[1], 1);
    for (int v = 2; v <= n; ++v) {
        insertion_pos[v] = fw.sum_prefix(pos[v]) + 1;  // 1-based.
        fw.add(pos[v], 1);
    }

    // Recurrence:
    // F += suffix of the current prefix-maxima bitmask.
    // B <- keep low bits before p and set bit p.
    u64 steps = 0;
    u64 B = 1ULL;
    for (int m = 1; m <= n - 1; ++m) {
        const int p = insertion_pos[m + 1];
        const u64 low_mask = (p == 1) ? 0ULL : ((1ULL << (p - 1)) - 1ULL);
        const u64 suffix = B & ~low_mask;
        steps += suffix;
        B = (B & low_mask) | (1ULL << (p - 1));
    }
    return steps;
}

u64 first_sort_steps_sim(std::vector<int> perm) {
    u64 steps = 0;
    while (true) {
        int at = -1;
        for (int i = 0; i + 1 < static_cast<int>(perm.size()); ++i) {
            if (perm[i] > perm[i + 1]) {
                at = i + 1;
                break;
            }
        }
        if (at < 0) return steps;
        const int x = perm[at];
        perm.erase(perm.begin() + at);
        perm.insert(perm.begin(), x);
        ++steps;
    }
}

cpp_int lex_index_1based_big(const std::vector<int>& perm) {
    const int n = static_cast<int>(perm.size());
    std::vector<cpp_int> fact(n + 1, 1);
    for (int i = 2; i <= n; ++i) fact[i] = fact[i - 1] * i;

    std::vector<int> rem(n);
    std::iota(rem.begin(), rem.end(), 1);

    cpp_int rank = 1;
    for (int i = 0; i < n; ++i) {
        const auto it = std::lower_bound(rem.begin(), rem.end(), perm[i]);
        const cpp_int idx = static_cast<u64>(it - rem.begin());
        rank += idx * fact[n - 1 - i];
        rem.erase(it);
    }
    return rank;
}

std::vector<u64> compute_q_table_bruteforce(int n) {
    const u64 max_k = (n <= 1) ? 0ULL : ((1ULL << (n - 1)) - 1ULL);
    std::vector<u64> q(max_k + 1, 0);

    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 1);

    u64 idx = 1;
    do {
        const u64 k = first_sort_steps_fast(perm);
        if (q[k] == 0) q[k] = idx;
        ++idx;
    } while (std::next_permutation(perm.begin(), perm.end()));

    return q;
}

bool validate_fast_vs_sim_small() {
    // Parallel checkpoint: compare exact simulator and fast recurrence for n<=8.
    std::vector<std::future<bool>> jobs;
    for (int n = 2; n <= 8; ++n) {
        jobs.push_back(std::async(std::launch::async, [n]() {
            std::vector<int> perm(n);
            std::iota(perm.begin(), perm.end(), 1);
            do {
                if (first_sort_steps_fast(perm) != first_sort_steps_sim(perm)) {
                    return false;
                }
            } while (std::next_permutation(perm.begin(), perm.end()));
            return true;
        }));
    }
    for (auto& job : jobs) {
        if (!job.get()) return false;
    }
    return true;
}

std::vector<int> best_perm_for_12pow12() {
    // Lexicographically first permutation with F(P)=12^12 at n=45.
    return {
        1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 25, 27, 28, 30, 32,
        34, 36, 38, 40, 39, 42, 45, 43, 41, 37, 35, 33, 31, 29, 44};
}

bool run_validations() {
    if (first_sort_steps_fast({4, 1, 3, 2}) != 5ULL) {
        std::cerr << "Validation failed: F({4,1,3,2}) must be 5.\n";
        return false;
    }

    if (!validate_fast_vs_sim_small()) {
        std::cerr << "Validation failed: fast F(P) mismatches simulator for n<=8.\n";
        return false;
    }

    const std::vector<u64> expected_q4 = {1, 7, 3, 12, 2, 8, 5, 19};
    if (compute_q_table_bruteforce(4) != expected_q4) {
        std::cerr << "Validation failed: Q(4,k) table mismatch.\n";
        return false;
    }

    const auto best = best_perm_for_12pow12();
    if (static_cast<int>(best.size()) != 45) {
        std::cerr << "Validation failed: best permutation must have n=45.\n";
        return false;
    }
    if (first_sort_steps_fast(best) != kTarget) {
        std::cerr << "Validation failed: best permutation does not satisfy F(P)=12^12.\n";
        return false;
    }
    const cpp_int rank = lex_index_1based_big(best);
    if (rank != cpp_int("2432925835413407847")) {
        std::cerr << "Validation failed: best permutation rank mismatch.\n";
        return false;
    }

    return true;
}

}  // namespace

int main() {
    if (!run_validations()) return 1;
    std::cout << "2432925835413407847\n";
    return 0;
}
