#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Segment {
    i64 x = 0;
    i64 y = 0;
    int dir = 0;
    u64 steps = 0;
};

struct Options {
    int order = 50;
    u64 steps = 1000000000000ULL;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--order=", options.order) ||
            parse_u64_after_prefix(arg, "--steps=", options.steps)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.order >= 0;
}

Segment combine(Segment lhs, Segment rhs) {
    for (int i = 0; i < lhs.dir; ++i) {
        const i64 old_x = rhs.x;
        rhs.x = rhs.y;
        rhs.y = -old_x;
    }

    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.dir = (lhs.dir + rhs.dir) & 3;
    lhs.steps += rhs.steps;
    return lhs;
}

struct Walker {
    std::vector<Segment> a;
    std::vector<Segment> b;
    Segment state;
    u64 remaining = 0;

    static constexpr Segment kF = {1, 0, 0, 1};
    static constexpr Segment kL = {0, 0, 3, 0};
    static constexpr Segment kR = {0, 0, 1, 0};

    explicit Walker(const int order, const u64 steps)
        : a(static_cast<std::size_t>(order + 1)), b(static_cast<std::size_t>(order + 1)), remaining(steps) {
        for (int i = order - 1; i >= 0; --i) {
            a[static_cast<std::size_t>(i)] = combine(
                combine(combine(combine(a[static_cast<std::size_t>(i + 1)], kR), b[static_cast<std::size_t>(i + 1)]), kF),
                kR);
            b[static_cast<std::size_t>(i)] = combine(
                combine(combine(combine(kL, kF), a[static_cast<std::size_t>(i + 1)]), kL),
                b[static_cast<std::size_t>(i + 1)]);
        }

        state = {0, 0, 3, 0};  // Facing North in this coordinate system.
    }

    void apply_segment(const Segment& seg) {
        state = combine(state, seg);
    }

    void execute_a(const int level) {
        if (remaining == 0 || level >= static_cast<int>(a.size()) - 1) {
            return;
        }

        const Segment& full = a[static_cast<std::size_t>(level)];
        if (remaining >= full.steps) {
            remaining -= full.steps;
            apply_segment(full);
            return;
        }

        execute_a(level + 1);
        if (remaining == 0) {
            return;
        }
        apply_segment(kR);

        execute_b(level + 1);
        if (remaining == 0) {
            return;
        }

        --remaining;
        apply_segment(kF);
        if (remaining == 0) {
            return;
        }

        apply_segment(kR);
    }

    void execute_b(const int level) {
        if (remaining == 0 || level >= static_cast<int>(b.size()) - 1) {
            return;
        }

        const Segment& full = b[static_cast<std::size_t>(level)];
        if (remaining >= full.steps) {
            remaining -= full.steps;
            apply_segment(full);
            return;
        }

        apply_segment(kL);

        if (remaining == 0) {
            return;
        }
        --remaining;
        apply_segment(kF);
        if (remaining == 0) {
            return;
        }

        execute_a(level + 1);
        if (remaining == 0) {
            return;
        }

        apply_segment(kL);
        execute_b(level + 1);
    }

    std::pair<i64, i64> run() {
        if (remaining > 0) {
            --remaining;
            apply_segment(kF);
        }
        execute_a(0);
        return {state.x, state.y};
    }
};

std::pair<i64, i64> solve(const int order, const u64 steps) {
    Walker walker(order, steps);
    return walker.run();
}

bool run_checkpoints() {
    const auto sample = solve(10, 500);
    if (sample.first != 18 || sample.second != 16) {
        std::cerr << "Checkpoint failed for D10 after 500 steps" << '\n';
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

    const auto answer = solve(options.order, options.steps);
    std::cout << answer.first << ',' << answer.second << '\n';
    return 0;
}
