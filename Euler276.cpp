#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    int perimeter_limit = 10'000'000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--perimeter=", options.perimeter_limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.perimeter_limit >= 3;
}

u64 all_triangles_with_perimeter(const int p) {
    if ((p & 1) == 0) {
        // round(p^2 / 48)
        return static_cast<u64>((static_cast<u64>(p) * static_cast<u64>(p) + 24ULL) / 48ULL);
    }

    // round((p + 3)^2 / 48)
    const u64 t = static_cast<u64>(p + 3);
    return (t * t + 24ULL) / 48ULL;
}

std::vector<int> mobius_up_to(const int n) {
    std::vector<int> mu(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> lp(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    mu[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (lp[static_cast<std::size_t>(i)] == 0) {
            lp[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(i) * static_cast<i64>(p);
            if (v > n || p > lp[static_cast<std::size_t>(i)]) {
                break;
            }
            lp[static_cast<std::size_t>(v)] = p;
            if (p == lp[static_cast<std::size_t>(i)]) {
                mu[static_cast<std::size_t>(v)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(v)] = -mu[static_cast<std::size_t>(i)];
        }
    }

    return mu;
}

u64 solve(const int perimeter_limit) {
    const std::vector<int> mu = mobius_up_to(perimeter_limit);

    std::vector<u64> prefix_all(static_cast<std::size_t>(perimeter_limit + 1), 0);
    for (int p = 1; p <= perimeter_limit; ++p) {
        prefix_all[static_cast<std::size_t>(p)] = prefix_all[static_cast<std::size_t>(p - 1)] +
                                                  all_triangles_with_perimeter(p);
    }

    __int128 sum = 0;
    for (int d = 1; d <= perimeter_limit; ++d) {
        const int mu_d = mu[static_cast<std::size_t>(d)];
        if (mu_d == 0) {
            continue;
        }
        const u64 add = prefix_all[static_cast<std::size_t>(perimeter_limit / d)];
        sum += static_cast<__int128>(mu_d) * static_cast<__int128>(add);
    }

    return static_cast<u64>(sum);
}

u64 brute_count(const int perimeter_limit) {
    u64 count = 0;
    for (int a = 1; a <= perimeter_limit / 3; ++a) {
        for (int b = a; b <= (perimeter_limit - a) / 2; ++b) {
            const int max_c = perimeter_limit - a - b;
            for (int c = b; c <= max_c; ++c) {
                if (a + b <= c) {
                    continue;
                }
                if (std::gcd(a, std::gcd(b, c)) != 1) {
                    continue;
                }
                ++count;
            }
        }
    }
    return count;
}

bool run_checkpoints() {
    for (int perimeter : {30, 50, 100, 150, 250}) {
        const u64 slow = brute_count(perimeter);
        const u64 fast = solve(perimeter);
        if (slow != fast) {
            std::cerr << "Checkpoint failed at perimeter=" << perimeter
                      << ": brute=" << slow << ", fast=" << fast << '\n';
            return false;
        }
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
    std::cout << solve(options.perimeter_limit) << '\n';
    return 0;
}
