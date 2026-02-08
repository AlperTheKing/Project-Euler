#include <array>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

struct Options {
    bool run_checkpoints = true;
};

using State = std::array<int, 4>;  // A2,A3,A4,A5 counts

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

int encode(const State& s) {
    return (((s[0] * 16 + s[1]) * 16 + s[2]) * 16 + s[3]);
}

double expected_single_envelopes(const State& state, std::unordered_map<int, double>& memo) {
    const int key = encode(state);
    auto it = memo.find(key);
    if (it != memo.end()) {
        return it->second;
    }

    const int total = state[0] + state[1] + state[2] + state[3];
    if (total == 0) {
        memo[key] = 0.0;
        return 0.0;
    }

    double ans = 0.0;

    // Count single-sheet states except the terminal single A5 case.
    if (total == 1 && state[3] == 0) {
        ans += 1.0;
    }

    for (int i = 0; i < 4; ++i) {
        if (state[static_cast<std::size_t>(i)] == 0) {
            continue;
        }

        State next = state;
        --next[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < 4; ++j) {
            ++next[static_cast<std::size_t>(j)];
        }

        const double probability = static_cast<double>(state[static_cast<std::size_t>(i)]) /
                                   static_cast<double>(total);
        ans += probability * expected_single_envelopes(next, memo);
    }

    memo[key] = ans;
    return ans;
}

double solve() {
    std::unordered_map<int, double> memo;
    const State start = {1, 1, 1, 1};
    return expected_single_envelopes(start, memo);
}

bool run_checkpoints() {
    std::unordered_map<int, double> memo;
    const State single_a4 = {0, 0, 1, 0};
    const double value = expected_single_envelopes(single_a4, memo);
    if (std::abs(value - 1.0) > 1e-12) {
        std::cerr << "Checkpoint failed for single A4 state" << '\n';
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

    std::cout << std::fixed << std::setprecision(6) << solve() << '\n';
    return 0;
}
