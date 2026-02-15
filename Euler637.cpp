#include <pthread.h>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;

struct BitMin {
    int n;
    std::vector<int> rightmost;
    std::vector<int>& arr;

    explicit BitMin(std::vector<int>& initial)
        : n(static_cast<int>(initial.size()) - 1), rightmost(8, -1), arr(initial) {}

    int suffix_min(int i) const {
        if (i == 0) return 0;
        for (int v = 0; v < static_cast<int>(rightmost.size()); ++v) {
            if (rightmost[static_cast<std::size_t>(v)] >= i) return v;
        }
        return n;
    }

    void update(int i, int new_val) {
        if (new_val >= arr[static_cast<std::size_t>(i)]) return;
        arr[static_cast<std::size_t>(i)] = new_val;
        seed(i, new_val);
    }

    void seed(int i, int val) {
        if (val < 0 || val >= n) return;
        if (val >= static_cast<int>(rightmost.size())) {
            int next_size = static_cast<int>(rightmost.size());
            while (next_size <= val) next_size <<= 1;
            rightmost.resize(static_cast<std::size_t>(next_size), -1);
        }
        int& pos = rightmost[static_cast<std::size_t>(val)];
        if (i > pos) pos = i;
    }
};

static std::vector<int> f_steps(int n, int b) {
    std::vector<int> steps(static_cast<std::size_t>(n + 1), 0);
    for (int i = b; i <= n; ++i) steps[static_cast<std::size_t>(i)] = n;
    std::vector<int> digit_sums(static_cast<std::size_t>(n + 1), 0);
    for (int i = 1; i <= n; ++i) {
        digit_sums[static_cast<std::size_t>(i)] = digit_sums[static_cast<std::size_t>(i / b)] + (i % b);
    }

    BitMin bit_min(steps);
    if (b > 1) bit_min.seed(b - 1, 0);

    std::vector<int> b_pows;
    for (i64 p = 1; p <= n; p *= b) b_pows.push_back(static_cast<int>(p));

    for (int k = b; k <= n; ++k) {
        const int dsum = digit_sums[static_cast<std::size_t>(k)];
        int current_min = steps[static_cast<std::size_t>(dsum)] + 1;

        auto dfs = [&](auto&& self, int krem, int digit_sum_rem, int exp, int curr_part, int part_sum) -> void {
            const int smallest = digit_sum_rem + part_sum + curr_part;
            const int cand = steps[static_cast<std::size_t>(smallest)] + 1;
            if (cand < current_min) current_min = cand;
            if (bit_min.suffix_min(smallest) + 1 >= current_min) return;
            if (krem == 0) return;

            const int digit = krem % b;
            const int next_k = krem / b;
            const int next_sum = digit_sum_rem - digit;
            self(self, next_k, next_sum, 1, digit, part_sum + curr_part);
            if (exp != 0 && exp < static_cast<int>(b_pows.size())) {
                self(self, next_k, next_sum, exp + 1,
                     curr_part + digit * b_pows[static_cast<std::size_t>(exp)], part_sum);
            }
        };
        dfs(dfs, k, dsum, 0, 0, 0);

        bit_min.update(k, current_min);
        steps[static_cast<std::size_t>(k)] = current_min;
    }

    return steps;
}

struct StepsTask {
    int n = 0;
    int b = 0;
    std::vector<int> steps;
};

static void* steps_worker(void* raw) {
    auto* task = static_cast<StepsTask*>(raw);
    task->steps = f_steps(task->n, task->b);
    return nullptr;
}

static i64 g(int n, int b1, int b2) {
    std::vector<int> s1;
    std::vector<int> s2;
    if (b1 == b2) {
        s1 = f_steps(n, b1);
        s2 = s1;
    } else {
        StepsTask t1{n, b1, {}};
        StepsTask t2{n, b2, {}};
        pthread_t th1{};
        pthread_t th2{};
        pthread_create(&th1, nullptr, steps_worker, &t1);
        pthread_create(&th2, nullptr, steps_worker, &t2);
        pthread_join(th1, nullptr);
        pthread_join(th2, nullptr);
        s1 = std::move(t1.steps);
        s2 = std::move(t2.steps);
    }

    i64 total = 0;
    for (int k = 1; k <= n; ++k) {
        if (s1[static_cast<std::size_t>(k)] == s2[static_cast<std::size_t>(k)]) total += k;
    }
    return total;
}

int main() {
    assert(g(100, 10, 3) == 3302);
    std::cout << g(10'000'000, 10, 3) << '\n';
    return 0;
}
