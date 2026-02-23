#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

namespace {

struct MatrixState {
    std::int64_t A;
    std::int64_t B;
    std::int64_t C;
    std::int64_t D;
};

struct Task {
    int m;
    int cap;
    MatrixState state;
    std::array<int, 12> prefix;
    int prefix_len;
    std::uint64_t local_count;
};

struct WorkerArgs {
    std::vector<Task>* tasks;
    std::atomic<std::size_t>* next_idx;
};

std::array<std::vector<int>, 13> g_proper_divisors;

std::array<std::vector<int>, 13> build_proper_divisors() {
    std::array<std::vector<int>, 13> divs;
    for (int m = 1; m <= 12; ++m) {
        for (int d = 1; d < m; ++d) {
            if (m % d == 0) {
                divs[m].push_back(d);
            }
        }
    }
    return divs;
}

MatrixState append_digit(const MatrixState& s, int digit) {
    const __int128 nextA = static_cast<__int128>(s.A) * digit + s.B;
    const __int128 nextB = -static_cast<__int128>(s.A);
    const __int128 nextC = static_cast<__int128>(s.C) * digit + s.D;
    const __int128 nextD = -static_cast<__int128>(s.C);
    return MatrixState{
        static_cast<std::int64_t>(nextA),
        static_cast<std::int64_t>(nextB),
        static_cast<std::int64_t>(nextC),
        static_cast<std::int64_t>(nextD)};
}

bool is_primitive(const std::array<int, 12>& seq, int m) {
    const std::vector<int>& divs = g_proper_divisors[m];
    for (int d : divs) {
        bool periodic = true;
        for (int i = d; i < m; ++i) {
            if (seq[i] != seq[i - d]) {
                periodic = false;
                break;
            }
        }
        if (periodic) {
            return false;
        }
    }
    return true;
}

void evaluate_m2(Task& task, std::array<int, 12>& seq) {
    const std::int64_t trace = task.state.A + task.state.D;
    if (trace < -1 || trace > 1) {
        return;
    }
    if (is_primitive(seq, 2)) {
        ++task.local_count;
    }
}

void evaluate_m3(Task& task, std::array<int, 12>& seq) {
    const MatrixState& s = task.state;
    const int cap = task.cap;
    std::array<unsigned char, 14> used{};
    for (int t = -1; t <= 1; ++t) {
        const __int128 rhs = static_cast<__int128>(t) - s.B + s.C;
        if (s.A == 0) {
            if (rhs != 0) {
                continue;
            }
            for (int x = 0; x <= cap; ++x) {
                if (used[static_cast<std::size_t>(x)] != 0) {
                    continue;
                }
                used[static_cast<std::size_t>(x)] = 1;
                seq[2] = x;
                if (is_primitive(seq, 3)) {
                    ++task.local_count;
                }
            }
            continue;
        }
        const __int128 den = s.A;
        if (rhs % den != 0) {
            continue;
        }
        const __int128 x128 = rhs / den;
        if (x128 < 0 || x128 > cap) {
            continue;
        }
        const int x = static_cast<int>(x128);
        if (used[static_cast<std::size_t>(x)] != 0) {
            continue;
        }
        used[static_cast<std::size_t>(x)] = 1;
        seq[2] = x;
        if (is_primitive(seq, 3)) {
            ++task.local_count;
        }
    }
}

void evaluate_leaf_suffix2(Task& task, const MatrixState& s, std::array<int, 12>& seq) {
    const int m = task.m;
    const int cap = task.cap;
    for (int u = 0; u <= cap; ++u) {
        seq[static_cast<std::size_t>(m - 2)] = u;
        const __int128 den = static_cast<__int128>(s.A) * u + s.B;
        const __int128 base = static_cast<__int128>(s.A) + s.D + static_cast<__int128>(s.C) * u;
        if (den == 0) {
            for (int t = -1; t <= 1; ++t) {
                if (base + t != 0) {
                    continue;
                }
                for (int v = 0; v <= cap; ++v) {
                    seq[static_cast<std::size_t>(m - 1)] = v;
                    if (is_primitive(seq, m)) {
                        ++task.local_count;
                    }
                }
            }
            continue;
        }
        for (int t = -1; t <= 1; ++t) {
            const __int128 num = base + t;
            if (num % den != 0) {
                continue;
            }
            const __int128 v128 = num / den;
            if (v128 < 0 || v128 > cap) {
                continue;
            }
            seq[static_cast<std::size_t>(m - 1)] = static_cast<int>(v128);
            if (is_primitive(seq, m)) {
                ++task.local_count;
            }
        }
    }
}

void dfs_prefix(Task& task, int depth, const MatrixState& s, std::array<int, 12>& seq) {
    if (depth == task.m - 2) {
        evaluate_leaf_suffix2(task, s, seq);
        return;
    }
    for (int digit = 0; digit <= task.cap; ++digit) {
        seq[static_cast<std::size_t>(depth)] = digit;
        const MatrixState next = append_digit(s, digit);
        dfs_prefix(task, depth + 1, next, seq);
    }
}

void process_task(Task& task) {
    std::array<int, 12> seq = task.prefix;
    task.local_count = 0;
    if (task.m == 2) {
        evaluate_m2(task, seq);
        return;
    }
    if (task.m == 3) {
        evaluate_m3(task, seq);
        return;
    }
    dfs_prefix(task, task.prefix_len, task.state, seq);
}

void* worker_entry(void* raw) {
    WorkerArgs* args = static_cast<WorkerArgs*>(raw);
    while (true) {
        const std::size_t idx = args->next_idx->fetch_add(1, std::memory_order_relaxed);
        if (idx >= args->tasks->size()) {
            break;
        }
        process_task((*args->tasks)[idx]);
    }
    return nullptr;
}

std::uint64_t count_exact_period(int m) {
    if (m == 1) {
        return 2;
    }

    const int cap = m + 1;
    const MatrixState id{1, 0, 0, 1};
    std::vector<Task> tasks;
    tasks.reserve(static_cast<std::size_t>(cap + 1) * static_cast<std::size_t>(cap + 1));
    for (int a0 = 0; a0 <= cap; ++a0) {
        const MatrixState s1 = append_digit(id, a0);
        for (int a1 = 0; a1 <= cap; ++a1) {
            const MatrixState s2 = append_digit(s1, a1);
            Task task{};
            task.m = m;
            task.cap = cap;
            task.state = s2;
            task.prefix.fill(0);
            task.prefix[0] = a0;
            task.prefix[1] = a1;
            task.prefix_len = 2;
            task.local_count = 0;
            tasks.push_back(task);
        }
    }

    long cpu = ::sysconf(_SC_NPROCESSORS_ONLN);
    std::size_t thread_count = 1;
    if (cpu > 0) {
        thread_count = static_cast<std::size_t>(cpu);
    }
    thread_count = std::max<std::size_t>(1, std::min<std::size_t>(thread_count, tasks.size()));

    std::atomic<std::size_t> next_idx{0};
    std::vector<pthread_t> threads(thread_count);
    std::vector<WorkerArgs> args(thread_count);
    for (std::size_t i = 0; i < thread_count; ++i) {
        args[i] = WorkerArgs{&tasks, &next_idx};
        pthread_create(&threads[i], nullptr, worker_entry, &args[i]);
    }
    for (std::size_t i = 0; i < thread_count; ++i) {
        pthread_join(threads[i], nullptr);
    }

    std::uint64_t total = 0;
    for (const Task& task : tasks) {
        total += task.local_count;
    }
    return total;
}

}  // namespace

int main() {
    g_proper_divisors = build_proper_divisors();

    std::array<std::uint64_t, 13> exact{};
    for (int m = 1; m <= 12; ++m) {
        exact[static_cast<std::size_t>(m)] = count_exact_period(m);
    }

    const std::uint64_t q1 = exact[1];
    const std::uint64_t q2 = exact[1] + exact[2];
    if (q1 != 2) {
        std::cerr << "Validation failed: Q(1) = " << q1 << ", expected 2\n";
        return 1;
    }
    if (q2 != 6) {
        std::cerr << "Validation failed: Q(2) = " << q2 << ", expected 6\n";
        return 1;
    }
    if (exact[2] != 4) {
        std::cerr << "Validation failed: exact period 2 = " << exact[2] << ", expected 4\n";
        return 1;
    }

    std::uint64_t q12 = 0;
    for (int m = 1; m <= 12; ++m) {
        q12 += exact[static_cast<std::size_t>(m)];
    }
    std::cout << q12 << '\n';
    return 0;
}
