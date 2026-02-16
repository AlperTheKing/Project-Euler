#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;

static int cycles_in_perm(const std::array<std::uint8_t, 24>& p) {
    bool seen[24] = {false};
    int cycles = 0;
    for (int i = 0; i < 24; ++i) {
        if (seen[i]) continue;
        ++cycles;
        int x = i;
        while (!seen[x]) {
            seen[x] = true;
            x = p[(size_t)x];
        }
    }
    return cycles;
}

static std::vector<std::array<std::uint8_t, 8>> build_orientations() {
    std::vector<std::array<std::uint8_t, 8>> oris;
    oris.reserve(2187);
    for (int mask = 0; mask < 2187; ++mask) {
        std::array<std::uint8_t, 8> ori{};
        int t = mask;
        int sum = 0;
        for (int i = 0; i < 7; ++i) {
            const int o = t % 3;
            t /= 3;
            ori[(size_t)i] = (std::uint8_t)o;
            sum += o;
        }
        ori[7] = (std::uint8_t)((3 - (sum % 3)) % 3);
        oris.push_back(ori);
    }
    return oris;
}

static std::vector<std::array<std::uint8_t, 8>> build_permutations() {
    std::vector<std::array<std::uint8_t, 8>> perms;
    perms.reserve(40320);
    std::array<std::uint8_t, 8> perm{};
    for (int i = 0; i < 8; ++i) perm[(size_t)i] = (std::uint8_t)i;
    do {
        perms.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));
    return perms;
}

struct HistWorkerCtx {
    const std::vector<std::array<std::uint8_t, 8>>* perms;
    const std::vector<std::array<std::uint8_t, 8>>* oris;
    std::size_t begin_idx;
    std::size_t end_idx;
    std::array<u64, 25> local_freq;
};

static void* hist_worker_main(void* ptr) {
    auto* ctx = static_cast<HistWorkerCtx*>(ptr);
    ctx->local_freq.fill(0);

    for (std::size_t pi = ctx->begin_idx; pi < ctx->end_idx; ++pi) {
        const auto& perm = (*ctx->perms)[pi];
        for (const auto& ori : *ctx->oris) {
            std::array<std::uint8_t, 24> p{};
            for (int pos = 0; pos < 8; ++pos) {
                const int c = perm[(size_t)pos];
                const int o = ori[(size_t)pos];
                for (int j = 0; j < 3; ++j) {
                    const int k = (j - o + 3) % 3;
                    p[(size_t)(3 * c + k)] = (std::uint8_t)(3 * pos + j);
                }
            }
            const int cyc = cycles_in_perm(p);
            ++ctx->local_freq[(size_t)cyc];
        }
    }

    return nullptr;
}

static unsigned choose_thread_count(std::size_t tasks) {
    long cpu_count = sysconf(_SC_NPROCESSORS_ONLN);
    unsigned threads = (cpu_count > 0) ? static_cast<unsigned>(cpu_count) : 1U;
    if (threads > 8U) {
        threads = 8U;
    }
    if (threads > tasks) {
        threads = static_cast<unsigned>(tasks);
    }
    if (threads == 0U) {
        threads = 1U;
    }
    return threads;
}

static std::array<u64, 25> cycle_histogram() {
    const auto oris = build_orientations();
    const auto perms = build_permutations();

    std::array<u64, 25> freq{};
    freq.fill(0);

    const unsigned thread_count = choose_thread_count(perms.size());
    if (thread_count == 1U) {
        HistWorkerCtx ctx{&perms, &oris, 0U, perms.size(), {}};
        hist_worker_main(&ctx);
        return ctx.local_freq;
    }

    std::vector<pthread_t> threads(thread_count);
    std::vector<HistWorkerCtx> ctx(thread_count);
    const std::size_t block = (perms.size() + thread_count - 1U) / thread_count;

    unsigned created = 0U;
    bool failed = false;
    for (unsigned t = 0; t < thread_count; ++t) {
        const std::size_t begin = static_cast<std::size_t>(t) * block;
        const std::size_t end = std::min(begin + block, perms.size());
        if (begin >= end) {
            ctx[t] = HistWorkerCtx{&perms, &oris, 0U, 0U, {}};
            continue;
        }
        ctx[t] = HistWorkerCtx{&perms, &oris, begin, end, {}};
        if (pthread_create(&threads[t], nullptr, hist_worker_main, &ctx[t]) != 0) {
            failed = true;
            break;
        }
        ++created;
    }

    for (unsigned t = 0; t < created; ++t) {
        pthread_join(threads[t], nullptr);
    }

    if (failed) {
        HistWorkerCtx seq{&perms, &oris, 0U, perms.size(), {}};
        hist_worker_main(&seq);
        return seq.local_freq;
    }

    for (unsigned t = 0; t < thread_count; ++t) {
        for (int k = 0; k <= 24; ++k) {
            freq[(size_t)k] += ctx[t].local_freq[(size_t)k];
        }
    }

    return freq;
}

static cpp_int orbits_from_hist(const std::array<u64, 25>& freq, int colours) {
    const u64 G = 88179840ULL;

    cpp_int sum = 0;
    cpp_int pw = 1;
    for (int k = 0; k <= 24; ++k) {
        if (freq[(size_t)k]) sum += cpp_int(freq[(size_t)k]) * pw;
        pw *= colours;
    }
    return sum / G;
}

int main() {
    const auto freq = cycle_histogram();
    if (orbits_from_hist(freq, 2) != 183) {
        std::cerr << "Validation failed for n=2\n";
        return 1;
    }

    std::cout << orbits_from_hist(freq, 10) << "\n";
    return 0;
}
