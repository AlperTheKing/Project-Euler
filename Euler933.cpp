#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <pthread.h>
#include <iostream>
#include <cstdint>
#include <unistd.h>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::int64_t;

constexpr u64 upper_bound_for_width(u64 width) {
    return 60 * width;
}

struct HeightTask {
    const std::vector<std::vector<std::uint8_t>>* grundy = nullptr;
    std::vector<std::uint8_t>* row = nullptr;
    std::vector<u64>* win_count = nullptr;
    u64 width = 0;
    u64 h_start = 0;
    u64 h_end = 0;
    u64 partial_sum = 0;
};

static void* compute_height_range(void* arg) {
    auto* task = static_cast<HeightTask*>(arg);
    const auto& g = *task->grundy;
    auto& out = *task->row;
    auto& counts = *task->win_count;
    const u64 w = task->width;
    u64 local_sum = 0;

    for (u64 h = task->h_start; h < task->h_end; ++h) {
        std::array<u32, 256> seen{};
        for (u64 a = 1; a <= (w - 1) / 2; ++a) {
            for (u64 b = 1; b <= (h - 1) / 2; ++b) {
                const u32 move_nim =
                    g[static_cast<size_t>(a)][static_cast<size_t>(b)] ^
                    g[static_cast<size_t>(w - a)][static_cast<size_t>(b)] ^
                    g[static_cast<size_t>(a)][static_cast<size_t>(h - b)] ^
                    g[static_cast<size_t>(w - a)][static_cast<size_t>(h - b)];
                ++seen[move_nim];
            }
        }

        if ((w & 1) == 0) {
            seen[0] += h - 1;
        }
        if ((h & 1) == 0) {
            seen[0] += w - 1;
        }
        if ((w & 1) == 0 && (h & 1) == 0) {
            --seen[0];
        }

        u32 mex = 0;
        while (mex < seen.size() && seen[mex] != 0) {
            ++mex;
        }
        assert(mex < seen.size());

        out[static_cast<size_t>(h)] = static_cast<std::uint8_t>(mex);
        counts[static_cast<size_t>(h - 2)] = seen[0];
        local_sum += seen[0];
    }

    task->partial_sum = local_sum;
    return nullptr;
}

u64 solve(u64 max_width, u64 max_height) {
    const u64 upper = upper_bound_for_width(max_width);
    std::vector<std::vector<std::uint8_t>> grundy(max_width + 1, std::vector<std::uint8_t>(upper + 1, 0));

    u64 total = 0;
    u64 n_bound = 1;
    u64 c_bound = 4;

    for (u64 width = 2; width <= max_width; ++width) {
        c_bound = std::max(c_bound, 2 * n_bound - 1);
        const u64 hh = std::min(c_bound, max_height);
        std::vector<u64> win_count(hh - 1);

        for (u64 h = 2; h <= hh; ++h) {
            win_count[h - 2] = h;
        }
        long cores = sysconf(_SC_NPROCESSORS_ONLN);
        if (cores < 1) cores = 1;
        const long task_count = std::min(cores, static_cast<long>(hh - 1));

        std::vector<HeightTask> tasks(static_cast<size_t>(task_count));
        std::vector<pthread_t> threads(static_cast<size_t>(task_count));
        const u64 per_thread = ((hh - 1) + static_cast<u64>(task_count) - 1) / static_cast<u64>(task_count);

        u64 total_wins = 0;
        for (long t = 0; t < task_count; ++t) {
            const u64 start = 2 + static_cast<u64>(t) * per_thread;
            const u64 end = std::min(hh + 1, start + per_thread);
            if (start >= end) {
                tasks[static_cast<size_t>(t)].partial_sum = 0;
                continue;
            }
            tasks[static_cast<size_t>(t)] = {
                &grundy,
                &grundy[static_cast<size_t>(width)],
                &win_count,
                width,
                start,
                end,
                0
            };
            pthread_create(&threads[static_cast<size_t>(t)], nullptr, compute_height_range, &tasks[static_cast<size_t>(t)]);
        }

        for (long t = 0; t < task_count; ++t) {
            if (threads[static_cast<size_t>(t)] != 0) {
                pthread_join(threads[static_cast<size_t>(t)], nullptr);
                total_wins += tasks[static_cast<size_t>(t)].partial_sum;
            }
        }

        total += total_wins;
        n_bound = hh;
        for (; n_bound > 0; --n_bound) {
            if (grundy[width][n_bound] != grundy[width][n_bound - 1]) {
                break;
            }
        }

        std::fill(grundy[width].begin() + hh + 1, grundy[width].end(), grundy[width][hh]);

        if (hh < max_height) {
            const u64 extra = max_height - hh;
            total += extra * (win_count.back() * 2 + (extra + 1) * (width - 1)) / 2;
        }
    }

    return total;
}

int main() {
    const auto start = std::chrono::high_resolution_clock::now();
    const u64 answer = solve(123, 1'234'567);
    const auto elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start);

    assert(answer > 0);
    std::cout << elapsed.count() << " s | " << answer << '\n';
    return 0;
}
