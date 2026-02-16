#include <cstdint>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int u = 100'000'000;
    int k = 100'000;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
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
        if (parse_int_after_prefix(arg, "--u=", options.u) ||
            parse_int_after_prefix(arg, "--k=", options.k)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.u >= 1 && options.k >= 1 && options.k <= options.u;
}

std::vector<std::uint16_t> divisor_counts(const int u) {
    std::vector<std::uint16_t> d(static_cast<std::size_t>(u + 1), 0U);
    std::vector<std::uint16_t> lp(static_cast<std::size_t>(u + 1), 0U);
    std::vector<std::uint8_t> exp(static_cast<std::size_t>(u + 1), 0U);
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(u / 10));

    d[1] = 1U;
    for (int i = 2; i <= u; ++i) {
        int minp = lp[static_cast<std::size_t>(i)];
        if (minp == 0) {
            primes.push_back(i);
            d[static_cast<std::size_t>(i)] = 2U;
            exp[static_cast<std::size_t>(i)] = 1U;
            minp = i;
        }

        for (int p : primes) {
            if (p > minp) {
                break;
            }
            const u64 x64 = static_cast<u64>(i) * static_cast<u64>(p);
            if (x64 > static_cast<u64>(u)) {
                break;
            }
            const int x = static_cast<int>(x64);

            lp[static_cast<std::size_t>(x)] = static_cast<std::uint16_t>(p);
            if (p == minp) {
                const std::uint8_t e = static_cast<std::uint8_t>(exp[static_cast<std::size_t>(i)] + 1U);
                exp[static_cast<std::size_t>(x)] = e;
                d[static_cast<std::size_t>(x)] = static_cast<std::uint16_t>(
                    d[static_cast<std::size_t>(i)] /
                    static_cast<std::uint16_t>(exp[static_cast<std::size_t>(i)] + 1U) *
                    static_cast<std::uint16_t>(e + 1U));
                break;
            }

            exp[static_cast<std::size_t>(x)] = 1U;
            d[static_cast<std::size_t>(x)] = static_cast<std::uint16_t>(d[static_cast<std::size_t>(i)] * 2U);
        }
    }
    return d;
}

u64 solve(const int u, const int k) {
    const std::vector<std::uint16_t> d = divisor_counts(u);
    std::deque<int> dq;
    u64 sum = 0ULL;

    for (int i = 1; i <= u; ++i) {
        while (!dq.empty() &&
               d[static_cast<std::size_t>(dq.back())] <= d[static_cast<std::size_t>(i)]) {
            dq.pop_back();
        }
        dq.push_back(i);
        while (!dq.empty() && dq.front() <= i - k) {
            dq.pop_front();
        }
        if (i >= k) {
            sum += static_cast<u64>(d[static_cast<std::size_t>(dq.front())]);
        }
    }
    return sum;
}

u64 brute(const int u, const int k) {
    const std::vector<std::uint16_t> d = divisor_counts(u);
    u64 ans = 0ULL;
    for (int n = 1; n <= u - k + 1; ++n) {
        std::uint16_t mx = 0U;
        for (int j = n; j <= n + k - 1; ++j) {
            if (d[static_cast<std::size_t>(j)] > mx) {
                mx = d[static_cast<std::size_t>(j)];
            }
        }
        ans += static_cast<u64>(mx);
    }
    return ans;
}

bool run_checkpoints() {
    if (solve(1000, 10) != 17'176ULL) {
        std::cerr << "Checkpoint failed: S(1000,10)=17176" << '\n';
        return false;
    }
    if (solve(2000, 37) != brute(2000, 37)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for (2000,37)" << '\n';
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
    std::cout << solve(options.u, options.k) << '\n';
    return 0;
}
