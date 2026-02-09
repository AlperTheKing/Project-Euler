#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr u64 kMod = 1'000'000'000ULL;

struct Options {
    u64 n = 450'283'905'890'997'363ULL;  // 3^37
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n == 0ULL) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

u64 normalize_mod(const i128 value) {
    i128 r = value % static_cast<i128>(kMod);
    if (r < 0) {
        r += static_cast<i128>(kMod);
    }
    return static_cast<u64>(r);
}

u64 reverse_bits_of_odd_core(u64 n) {
    while ((n & 1ULL) == 0ULL) {
        n >>= 1U;
    }
    u64 reversed = 0ULL;
    while (n > 0ULL) {
        reversed = (reversed << 1U) | (n & 1ULL);
        n >>= 1U;
    }
    return reversed;
}

class Solver {
public:
    Solver() {
        sum_cache_[0ULL] = 0ULL;
        sum_cache_[1ULL] = 1ULL;
        odd_cache_[0ULL] = 0ULL;
        odd_cache_[1ULL] = 1ULL;
    }

    u64 compute_sum_mod(u64 n) { return sum_mod(n); }

private:
    std::unordered_map<u64, u64> sum_cache_;
    std::unordered_map<u64, u64> odd_cache_;

    u64 f_mod(const u64 n) {
        return reverse_bits_of_odd_core(n) % kMod;
    }

    u64 sum_mod(const u64 n) {
        const auto it = sum_cache_.find(n);
        if (it != sum_cache_.end()) {
            return it->second;
        }

        u64 value = 0ULL;
        if ((n & 1ULL) == 0ULL) {
            const u64 half = n >> 1U;
            value = normalize_mod(static_cast<i128>(sum_mod(half)) + odd_sum_mod(half));
        } else {
            value = normalize_mod(static_cast<i128>(sum_mod(n - 1ULL)) + f_mod(n));
        }
        sum_cache_[n] = value;
        return value;
    }

    u64 odd_sum_mod(const u64 n) {
        const auto it = odd_cache_.find(n);
        if (it != odd_cache_.end()) {
            return it->second;
        }

        u64 value = 0ULL;
        if ((n & 1ULL) == 0ULL) {
            const u64 half = n >> 1U;
            value = normalize_mod(static_cast<i128>(-1) + 5 * static_cast<i128>(odd_sum_mod(half)) -
                                  3 * static_cast<i128>(sum_mod(half - 1ULL)));
        } else {
            const u64 half = n >> 1U;
            value = normalize_mod(static_cast<i128>(odd_sum_mod(n - 1ULL)) +
                                  2 * static_cast<i128>(f_mod(n)) -
                                  static_cast<i128>(f_mod(half)));
        }
        odd_cache_[n] = value;
        return value;
    }
};

bool run_checkpoints() {
    Solver solver;
    if (solver.compute_sum_mod(8ULL) != 22ULL) {
        std::cerr << "Checkpoint failed: S(8)\n";
        return false;
    }
    if (solver.compute_sum_mod(100ULL) != 3604ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }

    Solver solver;
    const u64 answer = solver.compute_sum_mod(options.n);
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}
