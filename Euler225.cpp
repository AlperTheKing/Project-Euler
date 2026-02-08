#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int target = 124;
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
        if (parse_int_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 1;
}

bool divides_tribonacci(const int n) {
    struct State {
        int a = 0;
        int b = 0;
        int c = 0;

        bool operator==(const State& other) const {
            return a == other.a && b == other.b && c == other.c;
        }
    };

    auto next = [&](const State& s) -> State {
        return {s.b, s.c, (s.a + s.b + s.c) % n};
    };

    const State start{1 % n, 1 % n, 1 % n};
    State tortoise = next(start);
    State hare = next(next(start));

    if (tortoise.c == 0 || hare.c == 0) {
        return true;
    }

    while (!(tortoise == hare)) {
        tortoise = next(tortoise);
        hare = next(next(hare));
        if (tortoise.c == 0 || hare.c == 0) {
            return true;
        }
    }

    State cur = tortoise;
    do {
        if (cur.c == 0) {
            return true;
        }
        cur = next(cur);
    } while (!(cur == tortoise));

    return false;
}

int solve(const int target) {
    int found = 0;
    int n = 1;

    while (found < target) {
        n += 2;
        if (!divides_tribonacci(n)) {
            ++found;
        }
    }

    return n;
}

bool run_checkpoints() {
    if (solve(1) != 27) {
        std::cerr << "Checkpoint failed for first odd non-divisor" << '\n';
        return false;
    }
    if (solve(10) != 221) {
        std::cerr << "Checkpoint failed for 10th odd non-divisor" << '\n';
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

    std::cout << solve(options.target) << '\n';
    return 0;
}
