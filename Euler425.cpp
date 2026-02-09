#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int limit = 10000000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2;
}

std::vector<uint8_t> sieve_is_prime(int n) {
    std::vector<uint8_t> is_prime(static_cast<std::size_t>(n + 1), 1);
    is_prime[0] = 0;
    is_prime[1] = 0;
    for (int p = 2; 1LL * p * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= n; q += p) {
            is_prime[static_cast<std::size_t>(q)] = 0;
        }
    }
    return is_prime;
}

i64 S_value(int limit) {
    const std::vector<uint8_t> is_prime = sieve_is_prime(limit);

    std::vector<int> primes;
    primes.reserve(limit / 10);
    for (int x = 2; x <= limit; ++x) {
        if (is_prime[static_cast<std::size_t>(x)]) {
            primes.push_back(x);
        }
    }

    std::vector<int> pow10{1};
    while (pow10.back() <= limit / 10) {
        pow10.push_back(pow10.back() * 10);
    }

    const int INF = std::numeric_limits<int>::max();
    std::vector<int> best(static_cast<std::size_t>(limit + 1), INF);
    if (limit < 2 || !is_prime[2]) {
        return 0;
    }

    using Node = std::pair<int, int>;  // (minimax cost, prime)
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    best[2] = 2;
    pq.push({2, 2});

    auto relax = [&](int from_cost, int v) {
        if (v < 2 || v > limit || !is_prime[static_cast<std::size_t>(v)]) {
            return;
        }
        const int cand = std::max(from_cost, v);
        if (cand < best[static_cast<std::size_t>(v)]) {
            best[static_cast<std::size_t>(v)] = cand;
            pq.push({cand, v});
        }
    };

    while (!pq.empty()) {
        const auto [cost, p] = pq.top();
        pq.pop();
        if (cost != best[static_cast<std::size_t>(p)]) {
            continue;
        }

        int len = 1;
        while (len < static_cast<int>(pow10.size()) && p >= pow10[len]) {
            ++len;
        }

        for (int pos = 0; pos < len; ++pos) {
            const int place = pow10[pos];
            const int old_digit = (p / place) % 10;
            for (int d = 0; d <= 9; ++d) {
                if (d == old_digit) {
                    continue;
                }
                if (pos == len - 1 && d == 0) {
                    continue;  // no leading zero
                }
                const int v = p + (d - old_digit) * place;
                relax(cost, v);
            }
        }

        if (len > 1) {
            const int v = p % pow10[len - 1];
            relax(cost, v);
        }

        if (len < static_cast<int>(pow10.size())) {
            const int base = pow10[len];
            for (int d = 1; d <= 9; ++d) {
                const i64 v = static_cast<i64>(d) * base + p;
                if (v > limit) {
                    break;
                }
                relax(cost, static_cast<int>(v));
            }
        }
    }

    i64 ans = 0;
    for (const int p : primes) {
        if (best[static_cast<std::size_t>(p)] > p) {
            ans += p;
        }
    }
    return ans;
}

bool run_checkpoints() {
    if (S_value(1000) != 431LL) {
        std::cerr << "Checkpoint failed: S(10^3)\n";
        return false;
    }
    if (S_value(10000) != 78728LL) {
        std::cerr << "Checkpoint failed: S(10^4)\n";
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

    std::cout << S_value(options.limit) << '\n';
    return 0;
}
