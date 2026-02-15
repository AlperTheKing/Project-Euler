#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <unistd.h>
#include <vector>

namespace {

constexpr int kK = 10;
constexpr int kN = 12;
constexpr std::uint64_t kMod = 1'234'567'891ULL;
constexpr std::uint64_t kLcgMul = 920'461ULL;
constexpr std::uint64_t kLcgAdd = 800'217'387'569ULL;
constexpr std::uint64_t kLcgMod = 1'000'000'000'000ULL;

struct Ranker {
    std::array<std::uint64_t, kN + 1> pow_k{};
    std::array<int, kN + 1> mu{};

    Ranker() {
        pow_k[0] = 1;
        for (int i = 1; i <= kN; ++i) {
            pow_k[i] = pow_k[i - 1] * static_cast<std::uint64_t>(kK);
        }
        for (int i = 1; i <= kN; ++i) {
            mu[i] = mobius(i);
        }
    }

    static int mobius(int x) {
        int n = x;
        int primes = 0;
        for (int p = 2; p * p <= n; ++p) {
            if (n % p != 0) {
                continue;
            }
            int exp = 0;
            while (n % p == 0) {
                n /= p;
                ++exp;
            }
            if (exp > 1) {
                return 0;
            }
            ++primes;
        }
        if (n > 1) {
            ++primes;
        }
        return (primes % 2 == 0) ? 1 : -1;
    }

    static int lyndon_prefix(int n, const std::array<int, kN + 1>& w) {
        int p = 1;
        for (int i = 2; i <= n; ++i) {
            if (w[i] < w[i - p]) {
                return p;
            }
            if (w[i] > w[i - p]) {
                p = i;
            }
        }
        return p;
    }

    static bool is_necklace(int n, const std::array<int, kN + 1>& w) {
        int p = 1;
        for (int i = 2; i <= n; ++i) {
            if (w[i] < w[i - p]) {
                return false;
            }
            if (w[i] > w[i - p]) {
                p = i;
            }
        }
        return (n % p == 0);
    }

    static void largest_necklace(int n,
                                 const std::array<int, kN + 1>& w,
                                 std::array<int, kN + 1>& neck) {
        neck = w;
        while (!is_necklace(n, neck)) {
            const int p = lyndon_prefix(n, neck);
            --neck[p];
            for (int i = p + 1; i <= n; ++i) {
                neck[i] = kK;
            }
        }
    }

    std::uint64_t T_value(int n, const std::array<int, kN + 1>& w) const {
        std::array<int, kN + 1> neck{};
        largest_necklace(n, w, neck);

        std::array<std::array<std::uint64_t, kN + 1>, kN + 1> B{};
        B[0][0] = 1;
        for (int t = 1; t <= n; ++t) {
            B[t][t] = 0;
            for (int j = t - 1; j >= 0; --j) {
                B[t][j] = B[t][j + 1] +
                          static_cast<std::uint64_t>(kK - neck[j + 1]) * B[t - j - 1][0];
            }
        }

        std::array<std::array<int, kN + 1>, kN + 1> suf{};
        for (int i = 2; i <= n; ++i) {
            int p = i - 1;
            for (int j = i; j <= n; ++j) {
                if (neck[j] > neck[j - p]) {
                    p = j;
                }
                suf[i][j] = j - p;
            }
        }

        std::uint64_t total = static_cast<std::uint64_t>(lyndon_prefix(n, neck));
        for (int t = 1; t <= n; ++t) {
            for (int j = 0; j < n; ++j) {
                if (j + t <= n) {
                    total += B[t - 1][0] * static_cast<std::uint64_t>(neck[j + 1] - 1) *
                             pow_k[n - t - j];
                } else {
                    int s = 0;
                    if (j >= n - t + 2) {
                        s = suf[n - t + 2][j];
                    }
                    if (neck[j + 1] > neck[s + 1]) {
                        total += B[n - j + s][s + 1] +
                                 static_cast<std::uint64_t>(neck[j + 1] - neck[s + 1] - 1) *
                                     B[n - j - 1][0];
                    }
                }
            }
        }
        return total;
    }

    std::uint64_t rank_lyndon(int n, const std::array<int, kN + 1>& w) const {
        std::int64_t r = 0;
        for (int i = 1; i <= n; ++i) {
            if (n % i != 0) {
                continue;
            }
            r += static_cast<std::int64_t>(mu[n / i]) *
                 static_cast<std::int64_t>(T_value(i, w));
        }
        return static_cast<std::uint64_t>(r / n);
    }

    std::uint64_t rank_db(int n, const std::array<int, kN + 1>& w) const {
        int t = 0;
        while (t + 1 <= n && w[t + 1] == kK) {
            ++t;
        }
        int j = t;
        while (j + 1 <= n && w[j + 1] == 1) {
            ++j;
        }
        if (t >= 1 && j == n) {
            return pow_k[n] - static_cast<std::uint64_t>(t) + 1ULL;
        }

        if (is_necklace(n, w)) {
            std::uint64_t r = 0;
            for (int i = 1; i <= n; ++i) {
                if (n % i == 0) {
                    r += static_cast<std::uint64_t>(i) * rank_lyndon(i, w);
                }
            }
            return 1ULL - static_cast<std::uint64_t>(lyndon_prefix(n, w)) + r;
        }

        std::array<int, kN + 1> neck = w;
        int s = 0;
        while (!is_necklace(n, neck)) {
            ++s;
            std::array<int, kN + 1> shifted{};
            for (int i = 1; i <= n; ++i) {
                const int idx = i + s;
                shifted[i] = (idx <= n) ? w[idx] : w[idx - n];
            }
            neck = shifted;
        }

        if (s != t) {
            return rank_db(n, neck) + static_cast<std::uint64_t>(lyndon_prefix(n, neck) - s);
        }
        if (lyndon_prefix(n, neck) < n) {
            return rank_db(n, neck) - static_cast<std::uint64_t>(s);
        }

        std::array<int, kN + 1> prev{};
        for (int i = n - s + 1; i <= n; ++i) {
            neck[i] = 1;
        }
        largest_necklace(n, neck, prev);
        return rank_db(n, prev) + static_cast<std::uint64_t>(lyndon_prefix(n, prev) - s);
    }
};

void to_digits1(std::uint64_t x, std::array<int, kN + 1>& w) {
    for (int i = kN; i >= 1; --i) {
        w[i] = static_cast<int>(x % 10ULL) + 1;
        x /= 10ULL;
    }
}

struct Item {
    std::uint64_t rank;
    std::uint64_t value;
};

int detect_thread_count() {
    const char* env = std::getenv("PE_THREADS");
    if (env != nullptr) {
        const int t = std::atoi(env);
        if (t > 0) {
            return t;
        }
    }
    const long nproc = sysconf(_SC_NPROCESSORS_ONLN);
    if (nproc > 0) {
        return static_cast<int>(nproc);
    }
    return 4;
}

struct RankTask {
    Item* items = nullptr;
    int begin = 0;
    int end = 0;
};

void* rank_worker(void* arg) {
    auto* task = static_cast<RankTask*>(arg);
    Ranker ranker;
    std::array<int, kN + 1> w{};
    for (int i = task->begin; i < task->end; ++i) {
        to_digits1(task->items[static_cast<std::size_t>(i)].value, w);
        task->items[static_cast<std::size_t>(i)].rank = ranker.rank_db(kN, w);
    }
    return nullptr;
}

std::uint64_t solve(int N, bool exact_small = false) {
    std::vector<Item> items(static_cast<std::size_t>(N));

    std::uint64_t a = 0;
    for (int i = 0; i < N; ++i) {
        a = (kLcgMul * a + kLcgAdd) % kLcgMod;
        items[static_cast<std::size_t>(i)].value = a;
    }

    int thread_count = detect_thread_count();
    if (thread_count > N) {
        thread_count = N;
    }
    if (thread_count < 1) {
        thread_count = 1;
    }

    std::vector<pthread_t> threads(static_cast<std::size_t>(thread_count));
    std::vector<RankTask> tasks(static_cast<std::size_t>(thread_count));
    int begin = 0;
    for (int t = 0; t < thread_count; ++t) {
        const int block = N / thread_count + (t < (N % thread_count) ? 1 : 0);
        const int end = begin + block;
        tasks[static_cast<std::size_t>(t)] = RankTask{items.data(), begin, end};
        pthread_create(&threads[static_cast<std::size_t>(t)], nullptr, rank_worker,
                       &tasks[static_cast<std::size_t>(t)]);
        begin = end;
    }
    for (int t = 0; t < thread_count; ++t) {
        pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
    }

    std::sort(items.begin(), items.end(), [](const Item& lhs, const Item& rhs) {
        return lhs.rank < rhs.rank;
    });

    if (exact_small) {
        unsigned __int128 exact = 0;
        for (int i = 0; i < N; ++i) {
            exact += static_cast<unsigned __int128>(i + 1) *
                     static_cast<unsigned __int128>(items[static_cast<std::size_t>(i)].value);
        }
        return static_cast<std::uint64_t>(exact);
    }

    std::uint64_t ans = 0;
    for (int i = 0; i < N; ++i) {
        const std::uint64_t pos = static_cast<std::uint64_t>(i + 1) % kMod;
        const std::uint64_t val = items[static_cast<std::size_t>(i)].value % kMod;
        ans = (ans + (pos * val) % kMod) % kMod;
    }
    return ans;
}

void run_validations() {
    assert(solve(2, true) == 2'194'210'461'325ULL);
    assert(solve(10, true) == 32'698'850'376'317ULL);
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve(10'000'000) << '\n';
    return 0;
}
