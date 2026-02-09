#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kLimit = 1000000000000000000ULL;
constexpr u64 kMod = 100000000ULL;

struct Precomputed {
    std::vector<u64> fib;        // fib[0]=1, fib[1]=1, fib[2]=2, ...
    std::vector<u64> prefix_len; // P_k
    std::vector<u128> tail_sum;  // tail sum in interval k
    std::vector<u128> interval_sum;
    std::vector<u128> interval_prefix_sum;
};

u128 tri(const u64 n) {
    return static_cast<u128>(n) * static_cast<u128>(n + 1ULL) / 2U;
}

u128 sum_range(const u64 a, const u64 b) {
    if (b < a) {
        return 0;
    }
    const u64 len = b - a + 1ULL;
    return static_cast<u128>(a + b) * static_cast<u128>(len) / 2U;
}

Precomputed build_precomputed(const u64 max_n) {
    Precomputed data;
    data.fib = {1ULL, 1ULL, 2ULL};  // index 0..2
    while (data.fib.back() <= max_n) {
        const std::size_t sz = data.fib.size();
        data.fib.push_back(data.fib[sz - 1] + data.fib[sz - 2]);
    }

    const std::size_t kmax = data.fib.size() - 2;  // largest interval index with fib[k] <= max_n

    data.prefix_len.assign(data.fib.size(), 0ULL);
    data.tail_sum.assign(data.fib.size(), 0U);
    data.interval_sum.assign(data.fib.size(), 0U);
    data.interval_prefix_sum.assign(data.fib.size(), 0U);

    // Base prefix lengths.
    data.prefix_len[1] = 1;
    data.prefix_len[2] = 1;
    data.prefix_len[3] = 2;
    data.prefix_len[4] = 3;
    data.prefix_len[5] = 4;

    for (std::size_t k = 6; k <= kmax; ++k) {
        data.prefix_len[k] = data.fib[k - 2] + data.prefix_len[k - 3];
    }

    // Tail sums for intervals.
    data.tail_sum[5] = 1;
    data.tail_sum[6] = 2;
    for (std::size_t k = 7; k <= kmax; ++k) {
        const u64 a = data.prefix_len[k - 3];
        const u64 b = data.prefix_len[k - 2] - 1ULL;
        data.tail_sum[k] = data.tail_sum[k - 2] + sum_range(a, b);
    }

    // Sum of M(n) over a full interval [F_k, F_{k+1}-1].
    data.interval_sum[1] = 0;
    data.interval_sum[2] = 0;
    for (std::size_t k = 3; k <= kmax; ++k) {
        const u64 p = data.prefix_len[k];
        data.interval_sum[k] = tri(p - 1ULL) + data.tail_sum[k];
    }

    for (std::size_t k = 1; k <= kmax; ++k) {
        data.interval_prefix_sum[k] = data.interval_prefix_sum[k - 1] + data.interval_sum[k];
    }

    return data;
}

u128 prefix_tail(const Precomputed& data, const std::size_t k, const u64 j) {
    if (j == static_cast<u64>(-1)) {
        return 0;
    }
    if (k <= 4) {
        return 0;
    }
    if (k == 5) {
        return 1;
    }
    if (k == 6) {
        return 2;
    }

    const u64 start = data.prefix_len[k - 3];
    const u64 end = data.prefix_len[k - 2] - 1ULL;
    const u64 block_len = end - start + 1ULL;

    if (j < block_len) {
        return sum_range(start, start + j);
    }
    return sum_range(start, end) + prefix_tail(data, k - 2, j - block_len);
}

u128 prefix_interval(const Precomputed& data, const std::size_t k, const u64 idx) {
    const u64 interval_len = data.fib[k - 1];
    if (idx >= interval_len - 1ULL) {
        return data.interval_sum[k];
    }
    if (k <= 2) {
        return 0;
    }

    const u64 p = data.prefix_len[k];
    if (idx < p) {
        return tri(idx);
    }

    const u64 tail_idx = idx - p;
    return tri(p - 1ULL) + prefix_tail(data, k, tail_idx);
}

u128 fast_sum_M(const u64 n, const Precomputed& data) {
    if (n == 0ULL) {
        return 0;
    }

    std::size_t k = 1;
    while (k + 1 < data.fib.size() && data.fib[k + 1] <= n) {
        ++k;
    }

    const u128 full_before = data.interval_prefix_sum[k - 1];
    const u64 idx = n - data.fib[k];
    return full_before + prefix_interval(data, k, idx);
}

u64 smallest_zeckendorf_term(u64 n, const std::vector<u64>& fib12) {
    std::size_t i = fib12.size() - 1;
    u64 smallest = 0;
    while (n > 0ULL) {
        while (fib12[i] > n) {
            --i;
        }
        n -= fib12[i];
        smallest = fib12[i];
        if (i >= 2) {
            i -= 2;
        } else {
            break;
        }
    }
    return smallest;
}

u64 slow_sum_M(const int limit) {
    std::vector<u64> fib12 = {1ULL, 2ULL};
    while (fib12.back() <= static_cast<u64>(limit)) {
        const std::size_t sz = fib12.size();
        fib12.push_back(fib12[sz - 1] + fib12[sz - 2]);
    }

    std::vector<u64> z(static_cast<std::size_t>(limit + 1), 0ULL);
    for (int i = 1; i <= limit; ++i) {
        z[static_cast<std::size_t>(i)] = smallest_zeckendorf_term(static_cast<u64>(i), fib12);
    }

    u64 total = 0;
    for (int n = 1; n <= limit; ++n) {
        u64 best = 0;
        for (int x = 1; x < n; ++x) {
            if (2ULL * static_cast<u64>(x) < z[static_cast<std::size_t>(n - x)]) {
                best = static_cast<u64>(x);
            }
        }
        total += best;
    }
    return total;
}

bool run_checkpoints() {
    const Precomputed data = build_precomputed(1000000ULL);

    const u64 slow_100 = slow_sum_M(100);
    if (slow_100 != 728ULL) {
        std::cerr << "Checkpoint failed: slow sum for n<=100\n";
        return false;
    }

    const u128 fast_100 = fast_sum_M(100ULL, data);
    if (fast_100 != 728ULL) {
        std::cerr << "Checkpoint failed: fast sum for n<=100\n";
        return false;
    }

    const u64 slow_1000 = slow_sum_M(1000);
    const u128 fast_1000 = fast_sum_M(1000ULL, data);
    if (fast_1000 != slow_1000) {
        std::cerr << "Checkpoint failed: slow/fast mismatch for n<=1000\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const Precomputed data = build_precomputed(kLimit);
    const u128 total = fast_sum_M(kLimit, data);
    const u64 answer = static_cast<u64>(total % static_cast<u128>(kMod));
    std::cout << answer << '\n';
    return 0;
}
