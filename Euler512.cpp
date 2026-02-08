#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
    u64 n = 500'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
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
    return options.n >= 1ULL;
}

class OddTotientSummatory {
public:
    u64 sum(const u64 n) {
        if (n == 0ULL) {
            return 0ULL;
        }
        auto it = memo_.find(n);
        if (it != memo_.end()) {
            return it->second;
        }

        i128 result = static_cast<i128>(odd_part_prefix_sum(n));
        for (u64 l = 2ULL; l <= n;) {
            const u64 q = n / l;
            const u64 r = n / q;
            result -= static_cast<i128>(r - l + 1ULL) * static_cast<i128>(sum(q));
            l = r + 1ULL;
        }

        const u64 out = static_cast<u64>(result);
        memo_.emplace(n, out);
        return out;
    }

private:
    static u64 odd_numbers_sum(const u64 x) {
        const u64 m = (x + 1ULL) / 2ULL;
        return m * m;
    }

    static u64 odd_part_prefix_sum(u64 n) {
        u64 total = 0ULL;
        while (n > 0ULL) {
            total += odd_numbers_sum(n);
            n >>= 1ULL;
        }
        return total;
    }

    std::unordered_map<u64, u64> memo_;
};

u64 solve(const u64 n) {
    OddTotientSummatory solver;
    return solver.sum(n);
}

u64 brute_phi_odd_sum(const u64 n) {
    auto phi = [](u64 x) {
        u64 result = x;
        u64 t = x;
        for (u64 p = 2ULL; p * p <= t; ++p) {
            if (t % p != 0ULL) {
                continue;
            }
            while (t % p == 0ULL) {
                t /= p;
            }
            result -= result / p;
        }
        if (t > 1ULL) {
            result -= result / t;
        }
        return result;
    };
    u64 s = 0ULL;
    for (u64 i = 1ULL; i <= n; i += 2ULL) {
        s += phi(i);
    }
    return s;
}

bool run_checkpoints() {
    if (solve(100ULL) != 2007ULL) {
        std::cerr << "Checkpoint failed: g(100)=2007" << '\n';
        return false;
    }
    if (solve(10'000ULL) != brute_phi_odd_sum(10'000ULL)) {
        std::cerr << "Checkpoint failed: brute-force odd-phi sum for n=10000" << '\n';
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
        return 2;
    }
    std::cout << solve(options.n) << '\n';
    return 0;
}
