#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

struct WidthState {
    u64 threshold = 1;
    std::uint16_t tail_grundy = 0;
    std::vector<std::uint16_t> prefix;  // valid for heights h < threshold, index by h
};

std::uint16_t value_at(const WidthState& st, u64 h) {
    if (h < st.threshold) {
        return st.prefix[static_cast<std::size_t>(h)];
    }
    return st.tail_grundy;
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct SolveResult {
    u128 D = 0;
    std::vector<std::vector<u64>> C_small;  // C_small[w][h] for h <= stored height, optional
};

struct HeightRangeTask {
    const std::vector<std::vector<std::uint16_t>>* ax = nullptr;
    std::vector<std::uint16_t>* g = nullptr;
    std::vector<u64>* c = nullptr;
    u64 w = 0;
    u64 h_begin = 0;
    u64 h_end = 0;
};

void* solve_height_range(void* raw) {
    auto* task = static_cast<HeightRangeTask*>(raw);
    auto& g = *task->g;
    auto& c = *task->c;
    const auto& ax = *task->ax;
    const u64 w = task->w;

    std::vector<int> seen(64, 0);
    int stamp = 1;

    for (u64 h = task->h_begin; h < task->h_end; ++h) {
        ++stamp;
        if (stamp == 0) {
            std::fill(seen.begin(), seen.end(), 0);
            stamp = 1;
        }

        u64 cw = 0;
        const u64 half = (h - 1) / 2;
        const bool even_h = (h & 1ULL) == 0ULL;

        for (u64 i = 1; i < w; ++i) {
            const auto& a = ax[i];
            for (u64 j = 1; j <= half; ++j) {
                const std::uint16_t x = static_cast<std::uint16_t>(a[j] ^ a[h - j]);
                if (x >= seen.size()) {
                    seen.resize(static_cast<std::size_t>(x) + 64, 0);
                }
                seen[x] = stamp;
                if (x == 0) {
                    cw += 2;
                }
            }
            if (even_h) {
                seen[0] = stamp;
                ++cw;
            }
        }

        std::uint16_t mex = 0;
        while (mex < seen.size() && seen[mex] == stamp) {
            ++mex;
        }
        g[h] = mex;
        c[h] = cw;
    }
    return nullptr;
}

SolveResult solve_game(u64 W, u64 H, bool keep_small_table) {
    std::vector<WidthState> states(W + 1);

    states[1].threshold = 1;
    states[1].tail_grundy = 0;
    states[1].prefix = {0};

    SolveResult out;
    if (keep_small_table) {
        out.C_small.assign(W + 1, std::vector<u64>(H + 1, 0));
    }

    for (u64 w = 2; w <= W; ++w) {
        u64 B_raw = 1;
        for (u64 i = 1; i < w; ++i) {
            const u64 tau = std::max(states[i].threshold, states[w - i].threshold);
            B_raw = std::max(B_raw, 2 * tau + 1);
        }

        const u64 limit = std::min(H, B_raw);

        std::vector<std::vector<std::uint16_t>> ax(w);
        for (u64 i = 1; i < w; ++i) {
            ax[i].assign(limit + 1, 0);
            for (u64 h = 1; h <= limit; ++h) {
                ax[i][h] = static_cast<std::uint16_t>(value_at(states[i], h) ^ value_at(states[w - i], h));
            }
        }

        std::vector<std::uint16_t> g(limit + 1, 0);
        std::vector<u64> c(limit + 1, 0);

        long cores = sysconf(_SC_NPROCESSORS_ONLN);
        if (cores < 1) {
            cores = 1;
        }
        const long task_count = std::min<long>(cores, static_cast<long>(limit));
        if (task_count <= 1) {
            HeightRangeTask task{&ax, &g, &c, w, 1, limit + 1};
            solve_height_range(&task);
        } else {
            std::vector<HeightRangeTask> tasks(static_cast<std::size_t>(task_count));
            std::vector<pthread_t> threads(static_cast<std::size_t>(task_count));
            std::vector<char> started(static_cast<std::size_t>(task_count), 0);
            const u64 per = (limit + static_cast<u64>(task_count) - 1) / static_cast<u64>(task_count);

            for (long t = 0; t < task_count; ++t) {
                const u64 begin = 1 + static_cast<u64>(t) * per;
                const u64 end = std::min(limit + 1, begin + per);
                if (begin >= end) {
                    continue;
                }
                tasks[static_cast<std::size_t>(t)] = HeightRangeTask{&ax, &g, &c, w, begin, end};
                pthread_create(&threads[static_cast<std::size_t>(t)], nullptr, solve_height_range,
                               &tasks[static_cast<std::size_t>(t)]);
                started[static_cast<std::size_t>(t)] = 1;
            }
            for (long t = 0; t < task_count; ++t) {
                if (started[static_cast<std::size_t>(t)]) {
                    pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
                }
            }
        }

        if (keep_small_table) {
            for (u64 h = 1; h <= limit; ++h) {
                if (h <= H) {
                    out.C_small[w][h] = c[h];
                }
            }
        }

        if (B_raw <= H) {
            std::vector<int> tail_seen(64, 0);
            const int tail_stamp = 1;
            tail_seen[0] = tail_stamp;

            i64 beta = 0;

            for (u64 i = 1; i < w; ++i) {
                const u64 tau = std::max(states[i].threshold, states[w - i].threshold);
                const std::uint16_t const_x = static_cast<std::uint16_t>(states[i].tail_grundy ^ states[w - i].tail_grundy);

                u64 z = 0;
                for (u64 t = 1; t < tau; ++t) {
                    const std::uint16_t edge = static_cast<std::uint16_t>(ax[i][t] ^ const_x);
                    if (edge >= tail_seen.size()) {
                        tail_seen.resize(static_cast<std::size_t>(edge) + 64, 0);
                    }
                    tail_seen[edge] = tail_stamp;
                    if (ax[i][t] == const_x) {
                        ++z;
                    }
                }

                beta += static_cast<i64>(-2 * static_cast<i64>(tau) + 1 + 2 * static_cast<i64>(z));
            }

            std::uint16_t g_inf = 0;
            while (g_inf < tail_seen.size() && tail_seen[g_inf] == tail_stamp) {
                ++g_inf;
            }

            assert(g[limit] == g_inf);
            const u64 c_tail_at_B = static_cast<u64>(static_cast<i64>(w - 1) * static_cast<i64>(B_raw) + beta);
            assert(c[limit] == c_tail_at_B);

            u64 T = limit;
            while (T >= 1 && g[T] == g_inf) {
                --T;
            }
            ++T;

            states[w].threshold = T;
            states[w].tail_grundy = g_inf;
            states[w].prefix.assign(T, 0);
            for (u64 h = 1; h < T; ++h) {
                states[w].prefix[h] = g[h];
            }

            for (u64 h = 2; h < B_raw; ++h) {
                out.D += c[h];
            }
            if (B_raw <= H) {
                const u128 cnt = static_cast<u128>(H - B_raw + 1);
                const u128 s = (static_cast<u128>(B_raw) + static_cast<u128>(H)) * cnt / 2;
                const i128 tail = static_cast<i128>(w - 1) * static_cast<i128>(s) +
                                  static_cast<i128>(beta) * static_cast<i128>(cnt);
                assert(tail >= 0);
                out.D += static_cast<u128>(tail);
            }
        } else {
            states[w].threshold = H + 1;
            states[w].tail_grundy = 0;
            states[w].prefix.assign(limit + 1, 0);
            for (u64 h = 1; h <= limit; ++h) {
                states[w].prefix[h] = g[h];
            }

            for (u64 h = 2; h <= H; ++h) {
                out.D += c[h];
            }
        }
    }

    return out;
}

void run_validations() {
    const SolveResult small = solve_game(12, 123, true);

    assert(small.C_small[5][3] == 4);
    assert(small.D == 327398);
}

}  // namespace

int main() {
    run_validations();
    const SolveResult ans = solve_game(123, 1'234'567, false);
    std::cout << to_string_u128(ans.D) << '\n';
    return 0;
}
