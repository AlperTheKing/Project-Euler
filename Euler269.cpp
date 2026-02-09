#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int max_power = 16;
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
        if (parse_int_after_prefix(arg, "--power=", options.max_power)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_power >= 1;
}

u64 pow10i(const int exp) {
    u64 value = 1;
    for (int i = 0; i < exp; ++i) {
        value *= 10;
    }
    return value;
}

u64 count_subset_for_length(const int length, const int last_digit, const std::vector<int>& roots) {
    if (roots.empty()) {
        return 0;
    }

    const int even_slots = (length + 1) / 2;
    const int msd_pos = length - 1;

    const std::vector<int> zero_state(static_cast<std::size_t>(roots.size()), 0);
    std::map<std::vector<int>, u64> states;
    states[zero_state] = 1;

    for (int j = 0; j < even_slots; ++j) {
        const int even_pos = 2 * j;
        const int odd_pos = 2 * j + 1;

        std::vector<int> even_digits;
        if (j == 0) {
            even_digits.push_back(last_digit);
        } else {
            even_digits.reserve(10);
            for (int d = 0; d <= 9; ++d) {
                even_digits.push_back(d);
            }
        }
        if (even_pos == msd_pos) {
            std::vector<int> filtered;
            filtered.reserve(even_digits.size());
            for (int d : even_digits) {
                if (d != 0) {
                    filtered.push_back(d);
                }
            }
            even_digits.swap(filtered);
        }

        std::vector<int> odd_digits;
        if (odd_pos >= length) {
            odd_digits.push_back(0);
        } else {
            odd_digits.reserve(10);
            for (int d = 0; d <= 9; ++d) {
                odd_digits.push_back(d);
            }
            if (odd_pos == msd_pos) {
                std::vector<int> filtered;
                filtered.reserve(odd_digits.size());
                for (int d : odd_digits) {
                    if (d != 0) {
                        filtered.push_back(d);
                    }
                }
                odd_digits.swap(filtered);
            }
        }

        std::map<std::vector<int>, u64> next_states;
        for (const auto& kv : states) {
            const std::vector<int>& carry = kv.first;
            const u64 ways = kv.second;

            for (int ed : even_digits) {
                for (int od : odd_digits) {
                    std::vector<int> next(carry.size(), 0);
                    bool ok = true;
                    for (std::size_t idx = 0; idx < roots.size(); ++idx) {
                        const int t = roots[idx];
                        const int denom = t * t;
                        const long long num = static_cast<long long>(t) * od + carry[idx] - ed;
                        if (num % denom != 0) {
                            ok = false;
                            break;
                        }
                        next[idx] = static_cast<int>(num / denom);
                    }
                    if (!ok) {
                        continue;
                    }
                    next_states[next] += ways;
                }
            }
        }

        states.swap(next_states);
        if (states.empty()) {
            return 0;
        }
    }

    const auto it = states.find(zero_state);
    return (it == states.end()) ? 0 : it->second;
}

u64 count_for_length(const int length) {
    u64 total = 0;

    // Root x=0: exactly numbers ending in zero.
    if (length >= 2) {
        total += 9ULL * pow10i(length - 2);
    }

    for (int d0 = 1; d0 <= 9; ++d0) {
        std::vector<int> divisors;
        for (int t = 1; t <= 9; ++t) {
            if (d0 % t == 0) {
                divisors.push_back(t);
            }
        }

        const int m = static_cast<int>(divisors.size());
        long long union_count = 0;
        for (int mask = 1; mask < (1 << m); ++mask) {
            std::vector<int> roots;
            roots.reserve(static_cast<std::size_t>(m));
            for (int i = 0; i < m; ++i) {
                if (((mask >> i) & 1) != 0) {
                    roots.push_back(divisors[static_cast<std::size_t>(i)]);
                }
            }

            const u64 cnt = count_subset_for_length(length, d0, roots);
            if ((__builtin_popcount(static_cast<unsigned int>(mask)) & 1) != 0) {
                union_count += static_cast<long long>(cnt);
            } else {
                union_count -= static_cast<long long>(cnt);
            }
        }

        total += static_cast<u64>(union_count);
    }

    return total;
}

u64 solve(const int max_power) {
    u64 total = 0;
    for (int len = 1; len <= max_power; ++len) {
        total += count_for_length(len);
    }

    // Add 10^max_power itself, which always has root x=0.
    total += 1;
    return total;
}

bool has_integer_root(u64 n) {
    std::vector<int> digits;
    while (n > 0) {
        digits.push_back(static_cast<int>(n % 10ULL));
        n /= 10ULL;
    }
    std::reverse(digits.begin(), digits.end());

    const int d0 = digits.back();
    std::vector<int> roots;
    if (d0 == 0) {
        roots.push_back(0);
    }
    for (int t = 1; t <= 9; ++t) {
        if (d0 % t == 0) {
            roots.push_back(-t);
        }
    }

    for (int r : roots) {
        long long value = 0;
        for (int d : digits) {
            value = value * static_cast<long long>(r) + d;
        }
        if (value == 0) {
            return true;
        }
    }
    return false;
}

u64 brute_count(const u64 limit) {
    u64 count = 0;
    for (u64 n = 1; n <= limit; ++n) {
        if (has_integer_root(n)) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(3) != brute_count(1000)) {
        std::cerr << "Checkpoint failed for limit 1000" << '\n';
        return false;
    }
    if (solve(4) != brute_count(10000)) {
        std::cerr << "Checkpoint failed for limit 10000" << '\n';
        return false;
    }
    if (solve(5) != 14696ULL) {
        std::cerr << "Checkpoint failed for Z(100000)" << '\n';
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
    std::cout << solve(options.max_power) << '\n';
    return 0;
}
