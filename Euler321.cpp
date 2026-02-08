#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
    int terms = 40;
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
        if (parse_int_after_prefix(arg, "--terms=", options.terms)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.terms >= 1;
}

struct PellState {
    i128 x;
    i128 y;
};

u64 sum_first_terms(const int count) {
    PellState a{5, 2};    // Gives n=1
    PellState b{11, 4};   // Gives n=3

    auto next_state = [](PellState s) {
        // Multiply x + y*sqrt(8) by 3 + sqrt(8).
        const i128 nx = 3 * s.x + 8 * s.y;
        const i128 ny = s.x + 3 * s.y;
        return PellState{nx, ny};
    };

    u64 sum = 0ULL;
    int taken = 0;
    while (taken < count) {
        const u64 na = static_cast<u64>(a.y - 1);
        const u64 nb = static_cast<u64>(b.y - 1);
        if (na < nb) {
            sum += na;
            a = next_state(a);
        } else {
            sum += nb;
            b = next_state(b);
        }
        ++taken;
    }
    return sum;
}

bool run_checkpoints() {
    if (sum_first_terms(5) != 99ULL) {
        std::cerr << "Checkpoint failed for first five-term sample sum" << '\n';
        return false;
    }
    if (sum_first_terms(1) != 1ULL) {
        std::cerr << "Checkpoint failed for first term" << '\n';
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
    std::cout << sum_first_terms(options.terms) << '\n';
    return 0;
}
