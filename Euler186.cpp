#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

struct Options {
    int subscribers = 1000000;
    int pm = 524287;
    int threshold_percent = 99;
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
        if (parse_int_after_prefix(arg, "--subscribers=", options.subscribers) ||
            parse_int_after_prefix(arg, "--pm=", options.pm) ||
            parse_int_after_prefix(arg, "--threshold-percent=", options.threshold_percent)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.subscribers >= 2 && options.pm >= 0 && options.pm < options.subscribers &&
           options.threshold_percent > 0 && options.threshold_percent <= 100;
}

class LaggedFibonacci {
public:
    explicit LaggedFibonacci(const int mod) : mod_(mod), idx_(0), generated_(0) {
        values_.assign(55, 0);
    }

    int next() {
        ++generated_;
        if (generated_ <= 55) {
            const std::int64_t k = static_cast<std::int64_t>(generated_);
            const std::int64_t value = 100003LL - 200003LL * k + 300007LL * k * k * k;
            const int result = static_cast<int>((value % mod_ + mod_) % mod_);
            values_[static_cast<std::size_t>(generated_ - 1)] = result;
            return result;
        }

        const int a = values_[idx_];
        const int b = values_[(idx_ + 31) % 55];  // k-24 in a 55-slot ring
        const int result = (a + b) % mod_;
        values_[idx_] = result;
        idx_ = (idx_ + 1) % 55;
        return result;
    }

private:
    int mod_;
    int idx_;
    int generated_;
    std::vector<int> values_;
};

struct DSU {
    std::vector<int> parent;
    std::vector<int> size;

    explicit DSU(const int n) : parent(static_cast<std::size_t>(n)), size(static_cast<std::size_t>(n), 1) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(const int x) {
        if (parent[static_cast<std::size_t>(x)] == x) {
            return x;
        }
        parent[static_cast<std::size_t>(x)] = find(parent[static_cast<std::size_t>(x)]);
        return parent[static_cast<std::size_t>(x)];
    }

    void unite(const int a, const int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) {
            return;
        }
        if (size[static_cast<std::size_t>(ra)] < size[static_cast<std::size_t>(rb)]) {
            std::swap(ra, rb);
        }
        parent[static_cast<std::size_t>(rb)] = ra;
        size[static_cast<std::size_t>(ra)] += size[static_cast<std::size_t>(rb)];
    }
};

int solve(const int subscribers, const int pm, const int threshold_percent) {
    LaggedFibonacci gen(subscribers);
    DSU dsu(subscribers);

    const int target_size = subscribers * threshold_percent / 100;
    int successful_calls = 0;
    while (dsu.size[static_cast<std::size_t>(dsu.find(pm))] < target_size) {
        const int caller = gen.next();
        const int called = gen.next();
        if (caller == called) {
            continue;
        }
        ++successful_calls;
        dsu.unite(caller, called);
    }
    return successful_calls;
}

bool run_checkpoints() {
    LaggedFibonacci gen(1000000);
    const int s1 = gen.next();
    const int s2 = gen.next();
    const int s3 = gen.next();
    const int s4 = gen.next();
    const int s5 = gen.next();
    const int s6 = gen.next();
    if (s1 != 200007 || s2 != 100053 || s3 != 600183 ||
        s4 != 500439 || s5 != 600863 || s6 != 701497) {
        std::cerr << "Checkpoint failed for Lagged Fibonacci prefix" << '\n';
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
    std::cout << solve(options.subscribers, options.pm, options.threshold_percent) << '\n';
    return 0;
}
