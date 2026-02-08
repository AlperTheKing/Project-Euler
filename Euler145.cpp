#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 limit = 1000000000ULL;
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

u64 pow10_u64(int e) {
    u64 v = 1;
    for (int i = 0; i < e; ++i) {
        v *= 10ULL;
    }
    return v;
}

int digits_count(u64 n) {
    int d = 0;
    do {
        ++d;
        n /= 10ULL;
    } while (n > 0);
    return d;
}

u64 count_for_length(const int length) {
    if (length <= 1) {
        return 0;
    }

    std::array<int, 19> outer_count{};
    std::array<int, 19> inner_count{};

    for (int a = 1; a <= 9; ++a) {
        for (int b = 1; b <= 9; ++b) {
            ++outer_count[static_cast<std::size_t>(a + b)];
        }
    }

    for (int a = 0; a <= 9; ++a) {
        for (int b = 0; b <= 9; ++b) {
            ++inner_count[static_cast<std::size_t>(a + b)];
        }
    }

    const int half = length / 2;
    std::array<int, 10> sums{};
    u64 total = 0;

    const auto apply_sequence = [&](const int center_twice, const u64 weight) -> u64 {
        int carry = 0;

        for (int i = 0; i < length; ++i) {
            int pair_sum = 0;
            if ((length % 2) == 1 && i == half) {
                pair_sum = center_twice;
            } else {
                const int idx = (i < half) ? i : (length - 1 - i);
                pair_sum = sums[static_cast<std::size_t>(idx)];
            }

            const int digit = (pair_sum + carry) % 10;
            if ((digit % 2) == 0) {
                return 0;
            }
            carry = (pair_sum + carry) / 10;
        }

        if (carry > 1) {
            return 0;
        }

        return weight;
    };

    const auto dfs = [&](auto&& self, int idx, u64 weight) -> void {
        if (idx == half) {
            if ((length % 2) == 0) {
                total += apply_sequence(0, weight);
            } else {
                for (int middle = 0; middle <= 9; ++middle) {
                    total += apply_sequence(2 * middle, weight);
                }
            }
            return;
        }

        for (int s = 0; s <= 18; ++s) {
            const int ways = (idx == 0) ? outer_count[static_cast<std::size_t>(s)]
                                        : inner_count[static_cast<std::size_t>(s)];
            if (ways == 0) {
                continue;
            }
            sums[static_cast<std::size_t>(idx)] = s;
            self(self, idx + 1, weight * static_cast<u64>(ways));
        }
    };

    dfs(dfs, 0, 1ULL);
    return total;
}

bool is_power_of_ten(const u64 x) {
    if (x == 0) {
        return false;
    }
    u64 v = x;
    while ((v % 10ULL) == 0ULL) {
        v /= 10ULL;
    }
    return v == 1ULL;
}

u64 solve(const u64 limit) {
    if (limit <= 1) {
        return 0;
    }

    const int max_len = digits_count(limit - 1);
    u64 total = 0;

    for (int len = 1; len < max_len; ++len) {
        total += count_for_length(len);
    }

    if (limit == pow10_u64(max_len)) {
        total += count_for_length(max_len);
        return total;
    }

    if (limit <= 10000000ULL) {
        for (u64 n = pow10_u64(max_len - 1); n < limit; ++n) {
            if ((n % 10ULL) == 0ULL) {
                continue;
            }
            u64 rev = 0;
            u64 t = n;
            while (t > 0) {
                rev = rev * 10ULL + (t % 10ULL);
                t /= 10ULL;
            }
            u64 s = n + rev;
            bool ok = true;
            while (s > 0) {
                if (((s % 10ULL) % 2ULL) == 0ULL) {
                    ok = false;
                    break;
                }
                s /= 10ULL;
            }
            if (ok) {
                ++total;
            }
        }
        return total;
    }

    throw std::runtime_error("Non power-of-10 limits above 10^7 are not supported");
}

bool run_checkpoints() {
    if (count_for_length(2) != 20ULL) {
        std::cerr << "Checkpoint failed for length 2" << '\n';
        return false;
    }
    if (count_for_length(3) != 100ULL) {
        std::cerr << "Checkpoint failed for length 3" << '\n';
        return false;
    }
    if (count_for_length(7) != 50000ULL) {
        std::cerr << "Checkpoint failed for length 7" << '\n';
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

    try {
        std::cout << solve(options.limit) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
