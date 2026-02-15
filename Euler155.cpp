#include <algorithm>
#include <climits>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int max_n = 18;
    bool run_checkpoints = true;
};

struct Fraction {
    u64 num = 0;
    u64 den = 1;

    bool operator<(const Fraction& other) const {
        if (num != other.num) {
            return num < other.num;
        }
        return den < other.den;
    }

    bool operator==(const Fraction& other) const {
        return num == other.num && den == other.den;
    }
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
        if (parse_int_after_prefix(arg, "--max-n=", options.max_n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.max_n >= 1;
}

Fraction make_fraction(u64 num, u64 den) {
    const u64 g = std::gcd(num, den);
    return {num / g, den / g};
}

Fraction parallel_combine(const Fraction& a, const Fraction& b) {
    const u128 num = static_cast<u128>(a.num) * b.den + static_cast<u128>(b.num) * a.den;
    const u128 den = static_cast<u128>(a.den) * b.den;
    return make_fraction(static_cast<u64>(num), static_cast<u64>(den));
}

u64 distinct_up_to(const int max_n) {
    std::vector<std::vector<Fraction>> ways(static_cast<std::size_t>(max_n + 1));
    std::vector<Fraction> all;
    all.reserve(4'000'000);

    ways[1] = {{1, 1}};
    all.push_back({1, 1});

    for (int n = 2; n <= max_n; ++n) {
        std::vector<Fraction> cur;
        std::size_t estimate = 0;
        for (int a = 1; a <= n / 2; ++a) {
            const int b = n - a;
            const auto& left = ways[static_cast<std::size_t>(a)];
            const auto& right = ways[static_cast<std::size_t>(b)];
            const std::size_t combos = left.size() * right.size();
            const std::size_t term_count = 2ULL * combos;
            if (estimate <= std::numeric_limits<std::size_t>::max() - term_count) {
                estimate += term_count;
            } else {
                estimate = std::numeric_limits<std::size_t>::max();
            }
        }
        cur.reserve(std::min<std::size_t>(estimate, 16'000'000));

        for (int a = 1; a <= n / 2; ++a) {
            const int b = n - a;
            const auto& left = ways[static_cast<std::size_t>(a)];
            const auto& right = ways[static_cast<std::size_t>(b)];

            for (std::size_t i = 0; i < left.size(); ++i) {
                const std::size_t j_start = (a == b) ? i : 0;
                for (std::size_t j = j_start; j < right.size(); ++j) {
                    const Fraction p = parallel_combine(left[i], right[j]);
                    cur.push_back(p);
                    if (p.num != p.den) {
                        cur.push_back({p.den, p.num});
                    }
                }
            }
        }

        std::sort(cur.begin(), cur.end());
        cur.erase(std::unique(cur.begin(), cur.end()), cur.end());

        ways[static_cast<std::size_t>(n)] = std::move(cur);
        all.insert(all.end(), ways[static_cast<std::size_t>(n)].begin(), ways[static_cast<std::size_t>(n)].end());
    }

    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    return static_cast<u64>(all.size());
}

bool run_checkpoints() {
    if (distinct_up_to(3) != 7ULL) {
        std::cerr << "Checkpoint failed for D(3)" << '\n';
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

    std::cout << distinct_up_to(options.max_n) << '\n';
    return 0;
}
