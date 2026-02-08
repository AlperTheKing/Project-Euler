#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

int next_state(const int x) {
    const int a = (x >> 5) & 1;
    const int b = (x >> 4) & 1;
    const int c = (x >> 3) & 1;
    const int g = a ^ (b & c);

    return ((x & 31) << 1) | g;
}

u64 count_cycle_assignments(const int len) {
    // Number of binary circular strings of length len with no adjacent ones:
    // F_{len-1} + F_{len+1}.
    std::vector<u64> fib(static_cast<std::size_t>(len + 2), 0ULL);
    fib[1] = 1ULL;
    for (int i = 2; i <= len + 1; ++i) {
        fib[static_cast<std::size_t>(i)] = fib[static_cast<std::size_t>(i - 1)] + fib[static_cast<std::size_t>(i - 2)];
    }
    return fib[static_cast<std::size_t>(len - 1)] + fib[static_cast<std::size_t>(len + 1)];
}

u64 solve() {
    std::array<std::uint8_t, 64> visited{};
    std::vector<int> cycle_lengths;

    for (int s = 0; s < 64; ++s) {
        if (visited[static_cast<std::size_t>(s)]) {
            continue;
        }
        int len = 0;
        int cur = s;
        while (!visited[static_cast<std::size_t>(cur)]) {
            visited[static_cast<std::size_t>(cur)] = 1;
            cur = next_state(cur);
            ++len;
        }
        cycle_lengths.push_back(len);
    }

    u64 answer = 1ULL;
    for (int len : cycle_lengths) {
        answer *= count_cycle_assignments(len);
    }
    return answer;
}

bool run_checkpoints() {
    if (count_cycle_assignments(3) != 4ULL) {
        std::cerr << "Checkpoint failed for cycle length 3" << '\n';
        return false;
    }

    std::array<std::uint8_t, 64> visited{};
    std::vector<int> lengths;
    for (int s = 0; s < 64; ++s) {
        if (visited[static_cast<std::size_t>(s)]) {
            continue;
        }
        int len = 0;
        int cur = s;
        while (!visited[static_cast<std::size_t>(cur)]) {
            visited[static_cast<std::size_t>(cur)] = 1;
            cur = next_state(cur);
            ++len;
        }
        lengths.push_back(len);
    }
    std::sort(lengths.begin(), lengths.end());
    const std::vector<int> expected = {1, 2, 3, 6, 6, 46};
    if (lengths != expected) {
        std::cerr << "Checkpoint failed for permutation cycle decomposition" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
