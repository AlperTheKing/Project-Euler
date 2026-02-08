#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <queue>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int target_left = 3;
    int target_below = 3;
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
        if (parse_int_after_prefix(arg, "--left=", options.target_left) ||
            parse_int_after_prefix(arg, "--below=", options.target_below)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.target_left >= 0 && options.target_below >= 0;
}

struct Node {
    int left = 0;
    int below = 0;
    long double x = 1.0L;
    long double y = 0.0L;
    long double side = 0.0L;
};

long double largest_side(const long double x, const long double y) {
    // Solve y + t = 1/(x + t) for t>0.
    return (std::sqrt((x - y) * (x - y) + 4.0L) - (x + y)) / 2.0L;
}

u64 binom(const int n, const int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    u64 numer = 1;
    u64 denom = 1;
    const int kk = std::min(k, n - k);
    for (int i = 1; i <= kk; ++i) {
        numer *= static_cast<u64>(n - kk + i);
        denom *= static_cast<u64>(i);
    }
    return numer / denom;
}

u64 solve(const int target_left, const int target_below) {
    const u64 target_occurrences = binom(target_left + target_below, target_left);

    auto cmp = [](const Node& a, const Node& b) {
        return a.side < b.side;
    };
    std::priority_queue<Node, std::vector<Node>, decltype(cmp)> pq(cmp);

    Node root;
    root.left = 0;
    root.below = 0;
    root.x = 1.0L;
    root.y = 0.0L;
    root.side = largest_side(root.x, root.y);
    pq.push(root);

    u64 index = 0;
    u64 seen = 0;
    u64 answer = 0;

    while (seen < target_occurrences) {
        const Node cur = pq.top();
        pq.pop();

        ++index;
        if (cur.left == target_left && cur.below == target_below) {
            ++seen;
            answer = index;
        }

        Node right;
        right.left = cur.left + 1;
        right.below = cur.below;
        right.x = cur.x + cur.side;
        right.y = cur.y;
        right.side = largest_side(right.x, right.y);
        pq.push(right);

        Node up;
        up.left = cur.left;
        up.below = cur.below + 1;
        up.x = cur.x;
        up.y = cur.y + cur.side;
        up.side = largest_side(up.x, up.y);
        pq.push(up);
    }

    return answer;
}

bool run_checkpoints() {
    if (solve(1, 1) != 50ULL) {
        std::cerr << "Checkpoint failed for index (1,1)" << '\n';
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

    std::cout << solve(options.target_left, options.target_below) << '\n';
    return 0;
}
