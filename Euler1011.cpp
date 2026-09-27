#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u32 TARGET = 3'000'000;

u32 modified_euclid(u32 a, u32 b) {
    while (a != 1 && b != 1) {
        if (a <= b) {
            b /= a;
        } else {
            a /= b;
        }
    }
    return a == 1 ? b : a;
}

u128 brute_force(const u32 n) {
    u128 total = 0;
    for (u32 a = 1; a < n; ++a) {
        for (u32 b = 1; b < n; ++b) total += modified_euclid(a, b);
    }
    return total;
}

struct Task {
    u32 n = 0;
    u32 root = 0;
    std::atomic<u32>* next = nullptr;
    u128 total = 0;
};

// Pairs a,q >= 2 with aq < n carry weight #{b < n : floor(b/a) = q} = min(a, n-aq);
// for m = min(a,q) and x >= m, f(x,m) = f(floor(x/m),m).
void* pair_worker(void* argument) {
    Task& task = *static_cast<Task*>(argument);
    const u32 n = task.n;
    std::vector<u32> value;
    for (u32 m = task.next->fetch_add(1); m <= task.root; m = task.next->fetch_add(1)) {
        const u32 limit = (n - 1) / m;
        if (value.size() <= limit) value.resize(limit + 1);
        value[1] = m;
        for (u32 x = 2; x < m; ++x) value[x] = modified_euclid(x, m);

        u128 sum = 0;
        for (u32 k = 1, x = m; x <= limit; ++k) {
            const u32 current = value[k];
            for (const u32 end = std::min(x + m - 1, limit); x <= end; ++x) {
                value[x] = current;
                const u32 rest = n - x * m;
                u64 weight = std::min(m, rest);
                if (x != m) weight += std::min(x, rest);
                sum += weight * current;
            }
        }
        task.total += sum;
    }
    return nullptr;
}

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

u128 solve(const u32 n, const unsigned thread_count) {
    const u128 triangle = static_cast<u128>(n) * (n - 1) / 2;
    u128 quotient_one = 0;
    for (u64 a = 2; a < n; ++a) quotient_one += a * (std::min<u64>(2 * a - 1, n - 1) - a);

    u32 root = 1;
    while (static_cast<u64>(root + 1) * (root + 1) < n) ++root;

    std::atomic<u32> next{2};
    std::vector<Task> tasks(thread_count);
    std::vector<pthread_t> threads(thread_count);
    for (unsigned t = 0; t < thread_count; ++t) {
        tasks[t].n = n;
        tasks[t].root = root;
        tasks[t].next = &next;
        require(pthread_create(&threads[t], nullptr, pair_worker, &tasks[t]) == 0, "pthread_create");
    }
    u128 paired = 0;
    for (unsigned t = 0; t < thread_count; ++t) {
        require(pthread_join(threads[t], nullptr) == 0, "pthread_join");
        paired += tasks[t].total;
    }
    return 3 * triangle - 2 + 2 * (quotient_one + paired);
}

std::string to_string(u128 value) {
    std::string digits;
    do {
        digits.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    } while (value != 0);
    std::reverse(digits.begin(), digits.end());
    return digits;
}

void run_tests(const unsigned thread_count) {
    require(modified_euclid(123, 456) == 3 && modified_euclid(456, 123) == 3, "f(123,456) = 3");
    require(brute_force(10) == 343 && solve(10, 1) == 343, "E(10) = 343");
    require(brute_force(100) == 269288 && solve(100, thread_count) == 269288, "E(100) = 269288");
    for (u32 n = 2; n <= 300; ++n) {
        require(solve(n, 1) == brute_force(n), "brute force n=" + std::to_string(n));
    }
    for (const u32 n : {1000U, 2023U, 4096U}) {
        require(solve(n, thread_count) == brute_force(n), "brute force n=" + std::to_string(n));
    }
    require(solve(200'003, 1) == solve(200'003, thread_count), "thread consistency");
    std::cout << "All checks passed.\n";
}

}  // namespace

int main() {
    try {
        const unsigned thread_count = std::min(16U, std::max(1U, std::thread::hardware_concurrency()));
        run_tests(thread_count);
        const auto start = std::chrono::steady_clock::now();
        const u128 answer = solve(TARGET, thread_count);
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        std::cout << to_string(answer) << '\n';
        std::cerr << "Computed in " << elapsed.count() << "s with " << thread_count << " threads.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
