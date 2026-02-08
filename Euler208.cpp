#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using State = std::array<i64, 5>;

struct Options {
    int steps = 70;
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
        if (parse_int_after_prefix(arg, "--steps=", options.steps)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.steps >= 0;
}

bool is_closed(const State& s) {
    return s[1] + s[4] == 0 &&
           s[2] + s[3] == 0 &&
           s[1] + s[2] == 0 &&
           s[0] + s[1] == 0;
}

u64 solve(const int steps) {
    std::map<State, u64> current;
    std::map<State, u64> next;
    current[{0, 0, 0, 0, 0}] = 1;

    for (int round = 0; round < steps; ++round) {
        next.clear();

        for (const auto& entry : current) {
            const State& s = entry.first;
            const u64 count = entry.second;

            const State clockwise = {s[2], s[3], s[4], -s[0], -s[1] - 1};
            const State anticlockwise = {-s[3], -s[4] + 1, s[0], s[1], s[2]};

            next[clockwise] += count;
            next[anticlockwise] += count;
        }

        current.swap(next);
    }

    u64 answer = 0;
    for (const auto& entry : current) {
        if (is_closed(entry.first)) {
            answer += entry.second;
        }
    }
    return answer;
}

bool run_checkpoints() {
    if (solve(25) != 70932ULL) {
        std::cerr << "Checkpoint failed for 25 steps" << '\n';
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

    std::cout << solve(options.steps) << '\n';
    return 0;
}
