#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using u128 = unsigned __int128;

namespace {

struct Term {
    int rank;
    u128 value;
};

struct Key {
    int steps_left;
    u128 x;
    bool prev_u;

    bool operator==(const Key& other) const {
        return steps_left == other.steps_left && x == other.x && prev_u == other.prev_u;
    }
};

struct KeyHash {
    std::size_t operator()(const Key& k) const {
        uint64_t lo = static_cast<uint64_t>(k.x);
        uint64_t hi = static_cast<uint64_t>(k.x >> 64);
        std::size_t h = std::hash<uint64_t>{}(lo) ^ (std::hash<uint64_t>{}(hi) << 1);
        h ^= std::hash<int>{}(k.steps_left + 0x9e3779b9);
        h ^= std::hash<int>{}(k.prev_u ? 1 : 0) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct PathCounter {
    std::unordered_map<Key, unsigned long long, KeyHash> memo;

    static bool is_pow2(u128 x) {
        return x != 0 && (x & (x - 1)) == 0;
    }

    // Count valid reverse paths from current value, ignoring ambiguity.
    unsigned long long count(int steps_left, u128 x, bool prev_u) {
        if (steps_left == 0) return 1;
        Key key{steps_left, x, prev_u};
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;

        unsigned long long total = 0;

        u128 y = x * 2;
        if (!is_pow2(y)) {
            total += count(steps_left - 1, y, false);
        }

        if (!prev_u && (x % 2 == 0) && (x % 3 == 1)) {
            y = (x - 1) / 3;
            if ((y % 2 == 1) && !is_pow2(y)) {
                total += count(steps_left - 1, y, true);
            }
        }

        memo.emplace(key, total);
        return total;
    }
};

struct Task {
    int steps_left;
    u128 x;
    bool prev_u;
    int d;
    int u;
    std::vector<Term> arr;
    bool amb;
};

std::vector<std::vector<int>> slope_rank;

int get_rank(int d, int u) {
    return slope_rank[d][u];
}

void init_slope_rank(int max_m) {
    using boost::multiprecision::cpp_int;
    std::vector<cpp_int> pow2(max_m + 1);
    std::vector<cpp_int> pow3(max_m + 1);
    pow2[0] = 1;
    pow3[0] = 1;
    for (int i = 1; i <= max_m; ++i) {
        pow2[i] = pow2[i - 1] * 2;
        pow3[i] = pow3[i - 1] * 3;
    }

    std::vector<std::pair<int, int>> pairs;
    pairs.reserve((max_m + 1) * max_m);
    for (int d = 0; d <= max_m; ++d) {
        for (int u = 1; u <= max_m; ++u) {
            pairs.emplace_back(d, u);
        }
    }

    auto cmp = [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        const cpp_int& a2 = pow2[a.first];
        const cpp_int& a3 = pow3[a.second];
        const cpp_int& b2 = pow2[b.first];
        const cpp_int& b3 = pow3[b.second];
        return a2 * b3 < b2 * a3;
    };
    std::sort(pairs.begin(), pairs.end(), cmp);

    slope_rank.assign(max_m + 1, std::vector<int>(max_m + 1, -1));
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        slope_rank[pairs[i].first][pairs[i].second] = i;
    }
}

bool is_pow2(u128 x) {
    return x != 0 && (x & (x - 1)) == 0;
}

unsigned long long dfs(int steps_left,
                       u128 x,
                       bool prev_u,
                       int d,
                       int u,
                       std::vector<Term>& arr,
                       bool amb,
                       PathCounter& counter) {
    if (amb) {
        return counter.count(steps_left, x, prev_u);
    }
    if (steps_left == 0) {
        return 0;
    }

    unsigned long long total = 0;

    auto handle = [&](u128 new_x, bool new_prev_u, int new_d, int new_u) {
        int new_rank = get_rank(new_d, new_u);
        int pos = 0;
        while (pos < static_cast<int>(arr.size()) && arr[pos].rank < new_rank) {
            ++pos;
        }

        bool new_amb = amb;
        if (!new_amb) {
            if (pos > 0 && new_x < arr[pos - 1].value) {
                new_amb = true;
            } else if (pos < static_cast<int>(arr.size()) && new_x > arr[pos].value) {
                new_amb = true;
            }
        }

        if (new_amb) {
            total += counter.count(steps_left - 1, new_x, new_prev_u);
            return;
        }

        arr.insert(arr.begin() + pos, Term{new_rank, new_x});
        total += dfs(steps_left - 1, new_x, new_prev_u, new_d, new_u, arr, new_amb, counter);
        arr.erase(arr.begin() + pos);
    };

    u128 y = x * 2;
    if (!is_pow2(y)) {
        handle(y, false, d + 1, u);
    }

    if (!prev_u && (x % 2 == 0) && (x % 3 == 1)) {
        y = (x - 1) / 3;
        if ((y % 2 == 1) && !is_pow2(y)) {
            handle(y, true, d, u + 1);
        }
    }

    return total;
}

void build_tasks(int depth,
                 int steps_left,
                 u128 x,
                 bool prev_u,
                 int d,
                 int u,
                 std::vector<Term>& arr,
                 bool amb,
                 std::vector<Task>& tasks,
                 unsigned long long& amb_sum,
                 PathCounter& counter) {
    if (amb) {
        amb_sum += counter.count(steps_left, x, prev_u);
        return;
    }
    if (depth == 0) {
        tasks.push_back(Task{steps_left, x, prev_u, d, u, arr, amb});
        return;
    }
    if (steps_left == 0) {
        return;
    }

    auto handle = [&](u128 new_x, bool new_prev_u, int new_d, int new_u) {
        int new_rank = get_rank(new_d, new_u);
        int pos = 0;
        while (pos < static_cast<int>(arr.size()) && arr[pos].rank < new_rank) {
            ++pos;
        }

        bool new_amb = amb;
        if (!new_amb) {
            if (pos > 0 && new_x < arr[pos - 1].value) {
                new_amb = true;
            } else if (pos < static_cast<int>(arr.size()) && new_x > arr[pos].value) {
                new_amb = true;
            }
        }

        if (new_amb) {
            amb_sum += counter.count(steps_left - 1, new_x, new_prev_u);
            return;
        }

        arr.insert(arr.begin() + pos, Term{new_rank, new_x});
        build_tasks(depth - 1, steps_left - 1, new_x, new_prev_u, new_d, new_u, arr, new_amb, tasks, amb_sum, counter);
        arr.erase(arr.begin() + pos);
    };

    u128 y = x * 2;
    if (!is_pow2(y)) {
        handle(y, false, d + 1, u);
    }

    if (!prev_u && (x % 2 == 0) && (x % 3 == 1)) {
        y = (x - 1) / 3;
        if ((y % 2 == 1) && !is_pow2(y)) {
            handle(y, true, d, u + 1);
        }
    }
}

unsigned long long compute_ambiguous(int m, unsigned threads) {
    if (m <= 1) return 0;

    std::vector<Term> arr;
    arr.reserve(m);
    // k = 2 gives a_m = 5 (smallest possible terminal value).
    arr.push_back(Term{get_rank(0, 1), static_cast<u128>(5)});

    std::vector<Task> tasks;
    tasks.reserve(4096);

    PathCounter counter;
    unsigned long long amb_sum = 0;

    int split_depth = 8;
    build_tasks(split_depth, m - 1, static_cast<u128>(5), false, 0, 1, arr, false, tasks, amb_sum, counter);

    if (tasks.empty()) return amb_sum;

    std::atomic<std::size_t> index(0);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    std::vector<unsigned long long> partial(threads, 0);

    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
                PathCounter local_counter;
                while (true) {
                    std::size_t i = index.fetch_add(1);
                    if (i >= tasks.size()) break;
                    const Task& task = tasks[i];
                    std::vector<Term> local_arr = task.arr;
                    partial[t] += dfs(task.steps_left, task.x, task.prev_u, task.d, task.u,
                                  local_arr, task.amb, local_counter);
                }
            });
        }
    for (auto& th : pool) th.join();

    for (auto v : partial) amb_sum += v;
    return amb_sum;
}

unsigned long long fibonacci(int n) {
    if (n <= 2) return 1;
    unsigned long long a = 1;
    unsigned long long b = 1;
    for (int i = 3; i <= n; ++i) {
        unsigned long long c = a + b;
        a = b;
        b = c;
    }
    return b;
}

} // namespace

int main() {
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;

    init_slope_rank(90);

    std::cout << "--- Validation Checkpoints ---\n";
    for (int m : {5, 10, 20}) {
        unsigned long long base = fibonacci(m);
        unsigned long long extra = compute_ambiguous(m, threads);
        unsigned long long f = base + extra;
        unsigned long long expected = (m == 5) ? 5ULL : (m == 10) ? 55ULL : 6771ULL;
        std::cout << "f(" << m << ") = " << f << (f == expected ? " [PASS]" : " [FAIL]") << "\n";
    }

    std::cout << "\n--- Final Solution ---\n";
    int target = 90;
    unsigned long long base = fibonacci(target);
    unsigned long long extra = compute_ambiguous(target, threads);
    unsigned long long result = base + extra;
    std::cout << "f(90) = " << result << "\n";

    return 0;
}
