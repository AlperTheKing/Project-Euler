#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 K = 1234567890ULL;
constexpr u64 OUTPUT_MOD = 1000000000000000000ULL;

u64 valuation_factorial(u64 n, const int p) {
    u64 sum = 0;
    const u64 prime = static_cast<u64>(p);
    while (n > 0) {
        n /= prime;
        sum += n;
    }
    return sum;
}

u64 minimal_n_for_target(const int p, const u64 target, const u64 previous) {
    if (target == 0ULL) {
        return 0ULL;
    }

    const u64 factor = static_cast<u64>(p - 1);
    const u64 lower = std::max(previous, static_cast<u64>(static_cast<u128>(target) * factor));

    u64 hi = lower;
    while (true) {
        const u64 have = valuation_factorial(hi, p);
        if (have >= target) {
            break;
        }
        hi = static_cast<u64>(static_cast<u128>(hi) + static_cast<u128>(target - have) * factor);
    }

    u64 lo = lower;
    while (lo < hi) {
        const u64 mid = lo + (hi - lo) / 2ULL;
        if (valuation_factorial(mid, p) >= target) {
            hi = mid;
        } else {
            lo = mid + 1ULL;
        }
    }
    return lo;
}

std::vector<int> build_spf(const int limit) {
    std::vector<int> spf(static_cast<std::size_t>(limit + 1));
    std::iota(spf.begin(), spf.end(), 0);

    for (int i = 2; static_cast<long long>(i) * i <= limit; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            if (spf[static_cast<std::size_t>(j)] == j) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return spf;
}

u128 solve_sum(const int upper) {
    const std::vector<int> spf = build_spf(upper);

    std::vector<u64> required(static_cast<std::size_t>(upper + 1), 0ULL);
    std::vector<u64> needed_n(static_cast<std::size_t>(upper + 1), 0ULL);

    u64 current_max = 0;
    u128 sum = 0;

    for (int i = 2; i <= upper; ++i) {
        int x = i;
        while (x > 1) {
            const int p = spf[static_cast<std::size_t>(x)];
            int exp = 0;
            while (x % p == 0) {
                x /= p;
                ++exp;
            }

            required[static_cast<std::size_t>(p)] += K * static_cast<u64>(exp);
            const u64 updated = minimal_n_for_target(p, required[static_cast<std::size_t>(p)],
                                                     needed_n[static_cast<std::size_t>(p)]);
            needed_n[static_cast<std::size_t>(p)] = updated;
            if (updated > current_max) {
                current_max = updated;
            }
        }

        if (i >= 10) {
            sum += static_cast<u128>(current_max);
        }
    }

    return sum;
}

bool run_checkpoints() {
    const u128 sample = solve_sum(1000);
    if (sample != static_cast<u128>(614538266565663ULL)) {
        std::cerr << "Checkpoint failed: S(1000) mismatch\n";
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

    const u128 total = solve_sum(1000000);
    const u64 answer = static_cast<u64>(total % static_cast<u128>(OUTPUT_MOD));
    std::cout << answer << '\n';
    return 0;
}
