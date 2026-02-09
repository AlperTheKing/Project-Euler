#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
constexpr int MOD = 1000000007;

struct Options {
    int n = 10000000;
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
    try {
        value = std::stoi(tail);
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
    return options.n >= 1;
}

int mod_pow(int base, i64 exp) {
    i64 result = 1;
    i64 cur = base % MOD;
    i64 e = exp;
    while (e > 0) {
        if (e & 1LL) {
            result = (result * cur) % MOD;
        }
        cur = (cur * cur) % MOD;
        e >>= 1LL;
    }
    return static_cast<int>(result);
}

std::vector<std::pair<int, int>> inadmissible_points(int n) {
    const int lim = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    const int lim2 = static_cast<int>(std::sqrt(static_cast<long double>(2LL * n)));

    std::vector<std::uint8_t> is_square(static_cast<std::size_t>(2 * n + 1), 0U);
    for (int c = 1; c <= lim2; ++c) {
        is_square[static_cast<std::size_t>(c * c)] = 1U;
    }

    std::vector<std::pair<int, int>> pts;
    pts.reserve(8000);
    for (int a = 1; a <= lim; ++a) {
        const int a2 = a * a;
        for (int b = 1; b <= lim; ++b) {
            const int b2 = b * b;
            if (is_square[static_cast<std::size_t>(a2 + b2)] != 0U) {
                pts.emplace_back(a2, b2);
            }
        }
    }

    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    return pts;
}

int P_mod(const int n) {
    const int maxv = 2 * n;
    std::vector<int> fact(static_cast<std::size_t>(maxv + 1), 1);
    for (int i = 1; i <= maxv; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<int>((1LL * fact[static_cast<std::size_t>(i - 1)] * i) % MOD);
    }

    std::vector<int> inv_fact(static_cast<std::size_t>(maxv + 1), 1);
    inv_fact[static_cast<std::size_t>(maxv)] = mod_pow(fact[static_cast<std::size_t>(maxv)], MOD - 2);
    for (int i = maxv; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<int>((1LL * inv_fact[static_cast<std::size_t>(i)] * i) % MOD);
    }

    auto nCr = [&](int nn, int rr) -> int {
        if (rr < 0 || rr > nn) {
            return 0;
        }
        return static_cast<int>(
            1LL * fact[static_cast<std::size_t>(nn)] * inv_fact[static_cast<std::size_t>(rr)] % MOD *
            inv_fact[static_cast<std::size_t>(nn - rr)] % MOD);
    };

    std::vector<std::pair<int, int>> pts = inadmissible_points(n);
    const int m = static_cast<int>(pts.size());
    std::vector<int> ways(static_cast<std::size_t>(m), 0);

    for (int i = 0; i < m; ++i) {
        const int xi = pts[static_cast<std::size_t>(i)].first;
        const int yi = pts[static_cast<std::size_t>(i)].second;
        i64 w = nCr(xi + yi, xi);

        for (int j = 0; j < i; ++j) {
            const int xj = pts[static_cast<std::size_t>(j)].first;
            const int yj = pts[static_cast<std::size_t>(j)].second;
            if (xj <= xi && yj <= yi) {
                const int paths = nCr((xi - xj) + (yi - yj), xi - xj);
                w -= 1LL * ways[static_cast<std::size_t>(j)] * paths % MOD;
                if (w < 0) {
                    w += MOD;
                }
            }
        }

        ways[static_cast<std::size_t>(i)] = static_cast<int>(w % MOD);
    }

    i64 ans = nCr(2 * n, n);
    for (int i = 0; i < m; ++i) {
        const int x = pts[static_cast<std::size_t>(i)].first;
        const int y = pts[static_cast<std::size_t>(i)].second;
        const int tail_paths = nCr((n - x) + (n - y), n - x);
        ans -= 1LL * ways[static_cast<std::size_t>(i)] * tail_paths % MOD;
        if (ans < 0) {
            ans += MOD;
        }
    }

    return static_cast<int>(ans % MOD);
}

bool run_checkpoints() {
    if (P_mod(5) != 252) {
        std::cerr << "Checkpoint failed: P(5)\n";
        return false;
    }
    if (P_mod(16) != 596994440) {
        std::cerr << "Checkpoint failed: P(16)\n";
        return false;
    }
    if (P_mod(1000) != 341920854) {
        std::cerr << "Checkpoint failed: P(1000)\n";
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

    std::cout << P_mod(options.n) << '\n';
    return 0;
}
