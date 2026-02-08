#include <array>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

constexpr int K = 2000;

struct Options {
    u64 index = 1000000000000000000ULL;
    int mod = 20092010;
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

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

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
        if (parse_u64_after_prefix(arg, "--index=", options.index) ||
            parse_int_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.mod >= 2;
}

using Poly = std::array<int, K>;

Poly multiply_mod(const Poly& a, const Poly& b, const int mod) {
    std::array<long long, 2 * K> tmp{};

    for (int i = 0; i < K; ++i) {
        const long long ai = a[static_cast<std::size_t>(i)];
        if (ai == 0) {
            continue;
        }
        for (int j = 0; j < K; ++j) {
            const int bj = b[static_cast<std::size_t>(j)];
            if (bj == 0) {
                continue;
            }
            tmp[static_cast<std::size_t>(i + j)] += ai * bj;
            tmp[static_cast<std::size_t>(i + j)] %= mod;
        }
    }

    for (int i = 2 * K - 1; i >= K; --i) {
        const long long v = tmp[static_cast<std::size_t>(i)] % mod;
        if (v == 0) {
            continue;
        }
        tmp[static_cast<std::size_t>(i - K)] += v;
        tmp[static_cast<std::size_t>(i - K + 1)] += v;
        tmp[static_cast<std::size_t>(i - K)] %= mod;
        tmp[static_cast<std::size_t>(i - K + 1)] %= mod;
    }

    Poly out{};
    for (int i = 0; i < K; ++i) {
        out[static_cast<std::size_t>(i)] = static_cast<int>(tmp[static_cast<std::size_t>(i)] % mod);
    }
    return out;
}

int solve(const u64 index, const int mod) {
    if (index < K) {
        return 1 % mod;
    }

    Poly acc{};
    Poly base{};
    acc[0] = 1;
    base[1] = 1;  // x

    u64 exp = index;
    while (exp > 0) {
        if (exp & 1ULL) {
            acc = multiply_mod(acc, base, mod);
        }
        exp >>= 1ULL;
        if (exp != 0) {
            base = multiply_mod(base, base, mod);
        }
    }

    int answer = 0;
    for (int c : acc) {
        answer += c;
        if (answer >= mod) {
            answer -= mod;
        }
    }
    return answer;
}

int solve_bruteforce(const int index, const int mod) {
    if (index < K) {
        return 1 % mod;
    }
    std::array<int, K> window{};
    window.fill(1 % mod);
    for (int i = K; i <= index; ++i) {
        const int next = (window[0] + window[1]) % mod;
        for (int j = 0; j < K - 1; ++j) {
            window[static_cast<std::size_t>(j)] = window[static_cast<std::size_t>(j + 1)];
        }
        window[K - 1] = next;
    }
    return window[K - 1];
}

bool run_checkpoints() {
    if (solve(0, 20092010) != 1 || solve(1999, 20092010) != 1) {
        std::cerr << "Checkpoint failed for initial terms" << '\n';
        return false;
    }
    if (solve(2000, 20092010) != 2 || solve(2001, 20092010) != 2) {
        std::cerr << "Checkpoint failed for first recursive terms" << '\n';
        return false;
    }

    const int test_index = 25000;
    const int fast = solve(test_index, 20092010);
    const int brute = solve_bruteforce(test_index, 20092010);
    if (fast != brute) {
        std::cerr << "Checkpoint failed for cross-check index 25000" << '\n';
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

    std::cout << solve(options.index, options.mod) << '\n';
    return 0;
}
