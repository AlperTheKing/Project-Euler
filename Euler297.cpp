#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    i64 limit = 100000000000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

class Solver {
public:
    Solver() {
        fib_.push_back(1);
        fib_.push_back(2);
        while (fib_.back() < 1000000000000000000LL) {
            const i64 n = fib_[fib_.size() - 1] + fib_[fib_.size() - 2];
            fib_.push_back(n);
        }
    }

    i64 solve(const i64 n) {
        memo_.clear();
        return S(n);
    }

private:
    i64 S(const i64 n) {
        if (n <= 4) {
            return n - 1;
        }
        const auto it = memo_.find(n);
        if (it != memo_.end()) {
            return it->second;
        }

        auto lb = std::lower_bound(fib_.begin(), fib_.end(), n);
        const i64 p = *(lb - 1);
        const i64 value = S(p) + (n - p) + S(n - p);
        memo_.emplace(n, value);
        return value;
    }

    std::vector<i64> fib_;
    std::unordered_map<i64, i64> memo_;
};

bool run_checkpoints() {
    Solver solver;
    if (solver.solve(1000000) != 7894453LL) {
        std::cerr << "Checkpoint failed for limit 1e6" << '\n';
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

    Solver solver;
    std::cout << solver.solve(options.limit) << '\n';
    return 0;
}
