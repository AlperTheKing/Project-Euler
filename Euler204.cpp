#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 limit = 1000000000ULL;
    int type = 100;
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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--type=", options.type)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1 && options.type >= 2;
}

std::vector<int> primes_up_to(const int n) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(n + 1), 1);
    is_prime[0] = 0;
    is_prime[1] = 0;

    for (int i = 2; i * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

u64 count_hamming(const u64 limit, const int type) {
    const std::vector<int> primes = primes_up_to(type);

    u64 count = 0;

    const auto dfs = [&](auto&& self, std::size_t idx, u64 value) -> void {
        if (idx == primes.size()) {
            ++count;
            return;
        }

        const u64 p = static_cast<u64>(primes[idx]);
        for (u64 x = value;;) {
            self(self, idx + 1, x);
            if (x > limit / p) {
                break;
            }
            x *= p;
        }
    };

    dfs(dfs, 0, 1ULL);
    return count;
}

bool run_checkpoints() {
    if (count_hamming(100000000ULL, 5) != 1105ULL) {
        std::cerr << "Checkpoint failed for type=5, limit=1e8" << '\n';
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

    std::cout << count_hamming(options.limit, options.type) << '\n';
    return 0;
}
