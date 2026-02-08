#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 120000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 3;
}

std::set<int> generate_distinct_sums(const int limit) {
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(limit + 1));

    const int m_max = static_cast<int>(std::sqrt(2.0L * limit)) + 4;
    for (int m = 2; m <= m_max; ++m) {
        for (int n = 1; n < m; ++n) {
            if (std::gcd(m, n) != 1) {
                continue;
            }
            if ((m - n) % 3 == 0) {
                continue;
            }

            const int x0 = m * m - n * n;
            const int y0 = n * (2 * m + n);
            if (x0 <= 0 || y0 <= 0) {
                continue;
            }

            if (x0 > limit || y0 > limit) {
                continue;
            }

            for (int k = 1; ; ++k) {
                const int x = k * x0;
                const int y = k * y0;
                if (x > limit || y > limit) {
                    break;
                }

                adj[static_cast<std::size_t>(x)].push_back(y);
                adj[static_cast<std::size_t>(y)].push_back(x);
            }
        }
    }

    for (auto& v : adj) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    std::set<int> sums;

    for (int p = 1; p <= limit; ++p) {
        const auto& vp = adj[static_cast<std::size_t>(p)];
        if (vp.size() < 2) {
            continue;
        }

        for (int q : vp) {
            if (q <= p) {
                continue;
            }
            if (p + q >= limit) {
                continue;
            }

            const auto& vq = adj[static_cast<std::size_t>(q)];

            auto it_p = std::upper_bound(vp.begin(), vp.end(), q);
            auto it_q = std::upper_bound(vq.begin(), vq.end(), q);

            while (it_p != vp.end() && it_q != vq.end()) {
                if (*it_p == *it_q) {
                    const int r = *it_p;
                    const int s = p + q + r;
                    if (s <= limit) {
                        sums.insert(s);
                    }
                    ++it_p;
                    ++it_q;
                } else if (*it_p < *it_q) {
                    ++it_p;
                } else {
                    ++it_q;
                }
            }
        }
    }

    return sums;
}

u64 solve(const int limit) {
    const std::set<int> sums = generate_distinct_sums(limit);
    u64 total = 0;
    for (int s : sums) {
        total += static_cast<u64>(s);
    }
    return total;
}

u64 brute_small(const int limit) {
    auto is_square = [](const int x) {
        const int r = static_cast<int>(std::sqrt(static_cast<long double>(x)));
        return r * r == x || (r + 1) * (r + 1) == x;
    };

    std::set<int> sums;
    for (int p = 1; p <= limit; ++p) {
        for (int q = p + 1; p + q <= limit; ++q) {
            if (!is_square(p * p + p * q + q * q)) {
                continue;
            }
            for (int r = q + 1; p + q + r <= limit; ++r) {
                if (is_square(p * p + p * r + r * r) && is_square(q * q + q * r + r * r)) {
                    sums.insert(p + q + r);
                }
            }
        }
    }

    u64 total = 0;
    for (int s : sums) {
        total += static_cast<u64>(s);
    }
    return total;
}

bool run_checkpoints() {
    const std::set<int> sums_1000 = generate_distinct_sums(1000);
    if (sums_1000.count(784) == 0) {
        std::cerr << "Checkpoint failed: 784 not generated for limit 1000" << '\n';
        return false;
    }
    if (solve(1000) != 784ULL) {
        std::cerr << "Checkpoint failed for limit 1000" << '\n';
        return false;
    }
    if (solve(1500) != brute_small(1500)) {
        std::cerr << "Checkpoint failed for brute-force cross-check limit 1500" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
