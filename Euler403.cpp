#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using i64 = long long;
using boost::multiprecision::cpp_int;

struct Options {
    i64 n = 1000000000000LL;
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
    try {
        value = std::stoll(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_i64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

cpp_int F1(const cpp_int& n) {
    return (n * n - n + 6) * (n + 1) / 6;
}

cpp_int F2(const cpp_int& n) {
    return (n * n - n + 12) * (n + 1) * (n + 2) / 24;
}

cpp_int F3(const cpp_int& n) {
    return (n * n - n + 20) * (n + 1) * (n + 2) * (n + 3) / 120;
}

cpp_int S_exact(i64 n) {
    const cpp_int N = n;
    cpp_int ans = 0;

    // Closed-form decomposition by discriminant class and floor-division blocks.
    ans += F2(N) + F2(N + 1) + F2(N - 2) + F1(N) + F1(N + 1);

    for (i64 x = 2; x < n;) {
        const i64 c = n / x;
        const i64 next = std::min(n / c + 1, n);

        ans += F3(cpp_int(next - 1 + c)) - F3(cpp_int(x - 1 + c));
        if (x > c) {
            ans -= F3(cpp_int(next - c - 2)) - F3(cpp_int(x - c - 2));
        } else {
            ans += F3(cpp_int(c - x)) - F3(cpp_int(c - next));
        }

        x = next;
    }

    return ans;
}

bool run_checkpoints() {
    if (S_exact(5) != 344) {
        std::cerr << "Checkpoint failed: S(5)\n";
        return false;
    }
    if (S_exact(100) != 26709528) {
        std::cerr << "Checkpoint failed: S(100)\n";
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

    cpp_int ans = S_exact(options.n) % cpp_int(100000000);
    if (ans < 0) {
        ans += cpp_int(100000000);
    }
    std::cout << ans << '\n';
    return 0;
}
