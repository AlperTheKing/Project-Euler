#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int n = 20'000'000;
    bool run_checkpoints = true;
};

class Fenwick {
public:
    explicit Fenwick(const int n) : bit_(static_cast<std::size_t>(n + 1), 0) {}

    void add(int idx, const int delta) {
        const int n = static_cast<int>(bit_.size()) - 1;
        while (idx <= n) {
            bit_[static_cast<std::size_t>(idx)] += delta;
            idx += idx & -idx;
        }
    }

    int sum(int idx) const {
        int result = 0;
        while (idx > 0) {
            result += bit_[static_cast<std::size_t>(idx)];
            idx -= idx & -idx;
        }
        return result;
    }

private:
    std::vector<int> bit_;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n <= 0) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

std::vector<std::int8_t> mobius_values(const int n) {
    std::vector<int> lp(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    std::vector<std::int8_t> mu(static_cast<std::size_t>(n + 1), 0);
    mu[1] = 1;

    for (int i = 2; i <= n; ++i) {
        if (lp[static_cast<std::size_t>(i)] == 0) {
            lp[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (const int p : primes) {
            const int64_t x = static_cast<int64_t>(p) * static_cast<int64_t>(i);
            if (x > n || p > lp[static_cast<std::size_t>(i)]) {
                break;
            }
            lp[static_cast<std::size_t>(x)] = p;
            if (p == lp[static_cast<std::size_t>(i)]) {
                mu[static_cast<std::size_t>(x)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(x)] =
                static_cast<std::int8_t>(-mu[static_cast<std::size_t>(i)]);
        }
    }
    return mu;
}

i64 count_positive_weighted_segments(
    const std::vector<std::int8_t>& mu,
    const int weight_pos,
    const int weight_neg) {
    const int n = static_cast<int>(mu.size()) - 1;
    std::vector<int> prefix(static_cast<std::size_t>(n + 1), 0);

    int cur = 0;
    int min_prefix = 0;
    int max_prefix = 0;
    for (int i = 1; i <= n; ++i) {
        const int m = static_cast<int>(mu[static_cast<std::size_t>(i)]);
        if (m == 1) cur += weight_pos;
        if (m == -1) cur += weight_neg;
        prefix[static_cast<std::size_t>(i)] = cur;
        if (cur < min_prefix) min_prefix = cur;
        if (cur > max_prefix) max_prefix = cur;
    }

    const i64 span = static_cast<i64>(max_prefix) - static_cast<i64>(min_prefix) + 1LL;
    constexpr i64 DENSE_LIMIT = 50'000'000LL;

    i64 positives = 0;
    if (span <= DENSE_LIMIT) {
        Fenwick bit(static_cast<int>(span) + 2);
        for (int i = 0; i <= n; ++i) {
            const int idx = prefix[static_cast<std::size_t>(i)] - min_prefix + 1;
            positives += static_cast<i64>(bit.sum(idx - 1));
            bit.add(idx, 1);
        }
        return positives;
    }

    std::vector<int> coords = prefix;
    std::sort(coords.begin(), coords.end());
    coords.erase(std::unique(coords.begin(), coords.end()), coords.end());

    Fenwick bit(static_cast<int>(coords.size()) + 2);
    for (int i = 0; i <= n; ++i) {
        const int v = prefix[static_cast<std::size_t>(i)];
        const int idx = static_cast<int>(
            std::lower_bound(coords.begin(), coords.end(), v) - coords.begin()) + 1;
        positives += static_cast<i64>(bit.sum(idx - 1));
        bit.add(idx, 1);
    }
    return positives;
}

i64 solve(const int n) {
    const std::vector<std::int8_t> mu = mobius_values(n);
    const i64 total_pairs = static_cast<i64>(n) * static_cast<i64>(n + 1) / 2LL;
    const i64 bad_a = count_positive_weighted_segments(mu, 99, -100);
    const i64 bad_b = count_positive_weighted_segments(mu, -100, 99);
    return total_pairs - bad_a - bad_b;
}

bool run_checkpoints() {
    if (solve(10) != 13LL) {
        std::cerr << "Checkpoint failed: C(10)\n";
        return false;
    }
    if (solve(500) != 16'676LL) {
        std::cerr << "Checkpoint failed: C(500)\n";
        return false;
    }
    if (solve(10'000) != 20'155'319LL) {
        std::cerr << "Checkpoint failed: C(10000)\n";
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
