#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using BigInt = boost::multiprecision::cpp_int;

struct Options {
    u64 p = 123456789ULL;
    u64 q = 987654321ULL;
    bool run_checkpoints = true;
};

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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--p=", options.p) ||
            parse_u64_after_prefix(arg, "--q=", options.q)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.p > 0 && options.q > 0;
}

std::vector<BigInt> build_f_values(const int limit) {
    std::vector<BigInt> f(static_cast<std::size_t>(limit + 1), 0);
    f[0] = 1;
    for (int i = 1; i <= limit; ++i) {
        if ((i & 1) == 1) {
            f[static_cast<std::size_t>(i)] = f[static_cast<std::size_t>(i / 2)];
        } else {
            f[static_cast<std::size_t>(i)] = f[static_cast<std::size_t>(i / 2)] +
                                             f[static_cast<std::size_t>(i / 2 - 1)];
        }
    }
    return f;
}

int smallest_n_for_ratio(const int p, const int q, const int search_limit) {
    const std::vector<BigInt> f = build_f_values(search_limit);
    for (int n = 1; n <= search_limit; ++n) {
        if (f[static_cast<std::size_t>(n)] * q == f[static_cast<std::size_t>(n - 1)] * p) {
            return n;
        }
    }
    return -1;
}

void shortened_expansion_dfs(const u64 p, const u64 q, std::vector<u64>& out) {
    if (q == 0) {
        return;
    }

    if (p > q) {
        const u64 next_p = (p - 1) % q + 1;
        shortened_expansion_dfs(next_p, q, out);
        out.push_back((p - next_p) / q);
    } else {
        const u64 next_q = q % p;
        shortened_expansion_dfs(p, next_q, out);
        out.push_back((q - next_q) / p);
    }
}

std::vector<u64> shortened_expansion(u64 p, u64 q) {
    const u64 g = std::gcd(p, q);
    p /= g;
    q /= g;

    std::vector<u64> out;
    shortened_expansion_dfs(p, q, out);
    return out;
}

std::vector<u64> shortened_binary_of(std::uint64_t n) {
    std::string bits;
    while (n > 0) {
        bits.push_back((n & 1ULL) != 0ULL ? '1' : '0');
        n >>= 1ULL;
    }
    std::reverse(bits.begin(), bits.end());

    std::vector<u64> out;
    if (bits.empty()) {
        return out;
    }

    u64 run = 1;
    for (std::size_t i = 1; i < bits.size(); ++i) {
        if (bits[i] == bits[i - 1]) {
            ++run;
        } else {
            out.push_back(run);
            run = 1;
        }
    }
    out.push_back(run);
    return out;
}

std::string as_csv(const std::vector<u64>& values) {
    std::string out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            out.push_back(',');
        }
        out += std::to_string(values[i]);
    }
    return out;
}

bool run_checkpoints() {
    const std::vector<BigInt> f = build_f_values(10);
    if (f[10] != 5) {
        std::cerr << "Checkpoint failed for f(10)" << '\n';
        return false;
    }

    const std::vector<u64> expected{4, 3, 1};
    const std::vector<u64> got = shortened_expansion(13, 17);
    if (got != expected) {
        std::cerr << "Checkpoint failed for shortened expansion of 13/17" << '\n';
        return false;
    }

    const int n = smallest_n_for_ratio(13, 17, 2000);
    if (n != 241) {
        std::cerr << "Checkpoint failed for smallest n with f(n)/f(n-1)=13/17" << '\n';
        return false;
    }
    if (shortened_binary_of(static_cast<u64>(n)) != expected) {
        std::cerr << "Checkpoint failed for shortened binary expansion of 241" << '\n';
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

    std::cout << as_csv(shortened_expansion(options.p, options.q)) << '\n';
    return 0;
}
