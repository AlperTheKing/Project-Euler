#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Options {
    int sides = 4;
    bool run_checkpoints = true;
};

struct Outcome {
    int pos = 0;
    double prob = 0.0;
    bool sent_to_jail = false;
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
        if (parse_int_after_prefix(arg, "--sides=", options.sides)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.sides >= 2;
}

int next_railway(const int pos) {
    if (pos < 5 || pos >= 35) {
        return 5;
    }
    if (pos < 15) {
        return 15;
    }
    if (pos < 25) {
        return 25;
    }
    return 35;
}

int next_utility(const int pos) {
    if (pos < 12 || pos >= 28) {
        return 12;
    }
    return 28;
}

std::vector<Outcome> resolve_square(const int pos) {
    if (pos == 30) {
        return {{10, 1.0, true}};
    }

    if (pos == 2 || pos == 17 || pos == 33) {
        std::vector<Outcome> out;
        out.push_back({0, 1.0 / 16.0, false});
        out.push_back({10, 1.0 / 16.0, true});
        out.push_back({pos, 14.0 / 16.0, false});
        return out;
    }

    if (pos == 7 || pos == 22 || pos == 36) {
        std::vector<Outcome> out;
        out.push_back({0, 1.0 / 16.0, false});
        out.push_back({10, 1.0 / 16.0, true});
        out.push_back({11, 1.0 / 16.0, false});
        out.push_back({24, 1.0 / 16.0, false});
        out.push_back({39, 1.0 / 16.0, false});
        out.push_back({5, 1.0 / 16.0, false});
        out.push_back({next_railway(pos), 1.0 / 16.0, false});
        out.push_back({next_railway(pos), 1.0 / 16.0, false});
        out.push_back({next_utility(pos), 1.0 / 16.0, false});

        // Go back 3 squares, then resolve again.
        const int back3 = (pos + 37) % 40;
        const auto back_out = resolve_square(back3);
        for (const auto& e : back_out) {
            out.push_back({e.pos, e.prob / 16.0, e.sent_to_jail});
        }

        out.push_back({pos, 6.0 / 16.0, false});
        return out;
    }

    return {{pos, 1.0, false}};
}

std::string solve(const int sides) {
    constexpr int kSquares = 40;
    constexpr int kDoubleStates = 3;
    constexpr int kStateCount = kSquares * kDoubleStates;

    std::array<std::array<double, kStateCount>, kStateCount> transition{};

    const double roll_prob = 1.0 / static_cast<double>(sides * sides);

    for (int pos = 0; pos < kSquares; ++pos) {
        for (int doubles = 0; doubles < kDoubleStates; ++doubles) {
            const int from = pos * kDoubleStates + doubles;

            for (int d1 = 1; d1 <= sides; ++d1) {
                for (int d2 = 1; d2 <= sides; ++d2) {
                    const bool rolled_double = (d1 == d2);
                    int next_doubles = rolled_double ? doubles + 1 : 0;

                    if (next_doubles >= 3) {
                        const int to = 10 * kDoubleStates;  // jail, doubles reset to 0
                        transition[from][to] += roll_prob;
                        continue;
                    }

                    const int moved = (pos + d1 + d2) % kSquares;
                    const auto outcomes = resolve_square(moved);
                    for (const auto& out : outcomes) {
                        const int final_doubles = out.sent_to_jail ? 0 : next_doubles;
                        const int to = out.pos * kDoubleStates + final_doubles;
                        transition[from][to] += roll_prob * out.prob;
                    }
                }
            }
        }
    }

    std::array<double, kStateCount> state{};
    state[0] = 1.0;

    for (int iter = 0; iter < 300; ++iter) {
        std::array<double, kStateCount> next{};
        for (int i = 0; i < kStateCount; ++i) {
            if (state[i] == 0.0) {
                continue;
            }
            for (int j = 0; j < kStateCount; ++j) {
                next[j] += state[i] * transition[i][j];
            }
        }
        state = next;
    }

    std::array<double, kSquares> square_prob{};
    for (int pos = 0; pos < kSquares; ++pos) {
        for (int d = 0; d < kDoubleStates; ++d) {
            square_prob[pos] += state[pos * kDoubleStates + d];
        }
    }

    std::vector<int> idx(kSquares);
    for (int i = 0; i < kSquares; ++i) {
        idx[static_cast<std::size_t>(i)] = i;
    }

    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        if (square_prob[a] != square_prob[b]) {
            return square_prob[a] > square_prob[b];
        }
        return a < b;
    });

    std::ostringstream out;
    out << std::setw(2) << std::setfill('0') << idx[0] << std::setw(2) << std::setfill('0') << idx[1]
        << std::setw(2) << std::setfill('0') << idx[2];
    return out.str();
}

bool run_checkpoints() {
    const std::string modal = solve(4);
    if (modal.size() != 6U) {
        std::cerr << "Checkpoint failed for modal string length" << '\n';
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

    std::cout << solve(options.sides) << '\n';
    return 0;
}
