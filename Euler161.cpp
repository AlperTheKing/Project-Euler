#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int width = 9;
    int height = 12;
    bool run_checkpoints = true;
};

int g_width = 0;
int g_state_count = 0;
std::vector<i64> g_cur;
std::vector<i64> g_next;

inline bool get_bit(const int mask, const int pos) {
    return ((mask >> pos) & 1) != 0;
}

inline int bit(const int pos) {
    return 1 << pos;
}

void dfs_transition(const int d, const int state, const int mask) {
    if (d < 0) {
        if (((mask + 1) & ((1 << g_width) - 1)) == 0) {
            g_next[static_cast<std::size_t>(mask >> g_width)] += g_cur[static_cast<std::size_t>(state)];
        }
        return;
    }

    dfs_transition(d - 1, state, mask);

    if (!get_bit(mask, d) && !get_bit(mask, d + g_width)) {
        dfs_transition(d - 1, state, mask | bit(d) | bit(d + g_width) | bit(d + 2 * g_width));
    }

    if (d >= 2) {
        dfs_transition(d - 3, state, mask | bit(d + 2 * g_width) | bit(d - 1 + 2 * g_width) | bit(d - 2 + 2 * g_width));
    }

    if (d >= 1 && !get_bit(mask, d + g_width) && !get_bit(mask, d + 2 * g_width - 1)) {
        dfs_transition(d - 2, state, mask | bit(d + g_width) | bit(d + 2 * g_width - 1) | bit(d + 2 * g_width));
    }

    if (d >= 1 && !get_bit(mask, d + g_width - 1) && !get_bit(mask, d + 2 * g_width - 1)) {
        dfs_transition(d - 2, state, mask | bit(d + g_width - 1) | bit(d + 2 * g_width - 1) | bit(d + 2 * g_width));
    }

    if (d >= 1 && !get_bit(mask, d + g_width - 1) && !get_bit(mask, d + g_width)) {
        dfs_transition(d - 1, state, mask | bit(d + g_width - 1) | bit(d + g_width) | bit(d + 2 * g_width));
    }

    if (d < g_width - 1 && !get_bit(mask, d + g_width) && !get_bit(mask, d + g_width + 1)) {
        dfs_transition(d - 1, state, mask | bit(d + g_width + 1) | bit(d + g_width) | bit(d + 2 * g_width));
    }
}

i64 solve(const int width, const int height) {
    if ((width * height) % 3 != 0) {
        return 0;
    }
    if (width <= 0 || width >= 11) {
        return 0;
    }

    g_width = width;
    g_state_count = 1 << (2 * width);
    g_cur.assign(static_cast<std::size_t>(g_state_count), 0);
    g_next.assign(static_cast<std::size_t>(g_state_count), 0);

    g_cur[static_cast<std::size_t>(g_state_count - 1)] = 1;

    for (int row = 0; row < height; ++row) {
        for (int s = 0; s < g_state_count; ++s) {
            if (g_cur[static_cast<std::size_t>(s)] == 0) {
                continue;
            }
            dfs_transition(width - 1, s, s);
        }
        g_cur.swap(g_next);
        std::fill(g_next.begin(), g_next.end(), 0);
    }

    return g_cur[static_cast<std::size_t>(g_state_count - 1)];
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
        if (parse_int_after_prefix(arg, "--width=", options.width) ||
            parse_int_after_prefix(arg, "--height=", options.height)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.width >= 1 && options.height >= 1;
}

bool run_checkpoints() {
    if (solve(2, 9) != 41LL) {
        std::cerr << "Checkpoint failed for 2x9" << '\n';
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

    std::cout << solve(options.width, options.height) << '\n';
    return 0;
}
