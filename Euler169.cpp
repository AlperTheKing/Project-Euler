#include <boost/multiprecision/cpp_int.hpp>

#include <cstdint>
#include <iostream>
#include <map>
#include <string>

namespace {

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;

struct Options {
    cpp_int n = cpp_int("10000000000000000000000000");
    bool run_checkpoints = true;
};

bool parse_cpp_int_after_prefix(const std::string& arg, const std::string& prefix, cpp_int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
    }

    value = cpp_int(tail);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_cpp_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.n >= 0;
}

u64 f(const cpp_int& n, std::map<cpp_int, u64>& memo) {
    auto it = memo.find(n);
    if (it != memo.end()) {
        return it->second;
    }

    u64 ans = 0;
    if (n == 0) {
        ans = 1;
    } else if ((n & 1) != 0) {
        ans = f((n - 1) >> 1, memo);
    } else {
        const cpp_int half = n >> 1;
        ans = f(half, memo) + f(half - 1, memo);
    }

    memo[n] = ans;
    return ans;
}

u64 solve(const cpp_int& n) {
    std::map<cpp_int, u64> memo;
    memo[cpp_int(0)] = 1;
    memo[cpp_int(-1)] = 0;
    return f(n, memo);
}

bool run_checkpoints() {
    if (solve(cpp_int(10)) != 5ULL) {
        std::cerr << "Checkpoint failed for f(10)" << '\n';
        return false;
    }
    if (solve(cpp_int(5)) != 2ULL) {
        std::cerr << "Checkpoint failed for f(5)" << '\n';
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
