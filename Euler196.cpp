#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = long long;

struct Options {
    int row_a = 5678027;
    int row_b = 7208785;
    int single_row = 0;
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
        if (parse_int_after_prefix(arg, "--row-a=", options.row_a) ||
            parse_int_after_prefix(arg, "--row-b=", options.row_b) ||
            parse_int_after_prefix(arg, "--row=", options.single_row)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.single_row != 0) {
        return options.single_row >= 2;
    }
    return options.row_a >= 2 && options.row_b >= 2;
}

u64 row_start(const int r) {
    return static_cast<u64>(r) * static_cast<u64>(r - 1) / 2ULL + 1ULL;
}

std::vector<int> small_primes_up_to(const int n) {
    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;

    for (int i = 2; i <= n; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
        for (int p : primes) {
            const long long v = 1LL * i * p;
            if (v > n) {
                break;
            }
            is_composite[static_cast<std::size_t>(v)] = 1;
            if (i % p == 0) {
                break;
            }
        }
    }

    return primes;
}

u64 solve_row(const int n) {
    const int r_min = n - 2;
    const int r_max = n + 2;

    const u64 start_value = row_start(r_min);
    const u64 end_value = row_start(r_max + 1) - 1;

    const u64 range_len = end_value - start_value + 1;
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(range_len), 1);

    if (start_value == 1ULL) {
        is_prime[0] = 0;
    }

    const int max_p = static_cast<int>(std::sqrt(static_cast<long double>(end_value))) + 1;
    const std::vector<int> primes = small_primes_up_to(max_p);

    for (int p : primes) {
        const u64 pp = static_cast<u64>(p) * static_cast<u64>(p);
        if (pp > end_value) {
            break;
        }

        u64 first = (start_value + static_cast<u64>(p) - 1ULL) / static_cast<u64>(p);
        first *= static_cast<u64>(p);
        if (first < pp) {
            first = pp;
        }

        for (u64 x = first; x <= end_value; x += static_cast<u64>(p)) {
            is_prime[static_cast<std::size_t>(x - start_value)] = 0;
        }
    }

    std::vector<std::vector<std::uint8_t>> row_prime(5);
    std::vector<std::vector<std::uint8_t>> row_good(5);

    for (int idx = 0; idx < 5; ++idx) {
        const int row = r_min + idx;
        row_prime[static_cast<std::size_t>(idx)].assign(static_cast<std::size_t>(row), 0);
        row_good[static_cast<std::size_t>(idx)].assign(static_cast<std::size_t>(row), 0);

        const u64 rs = row_start(row);
        for (int col = 0; col < row; ++col) {
            const u64 value = rs + static_cast<u64>(col);
            row_prime[static_cast<std::size_t>(idx)][static_cast<std::size_t>(col)] =
                is_prime[static_cast<std::size_t>(value - start_value)];
        }
    }

    for (int rid = 1; rid <= 3; ++rid) {
        const int row_len = r_min + rid;
        for (int col = 0; col < row_len; ++col) {
            if (!row_prime[static_cast<std::size_t>(rid)][static_cast<std::size_t>(col)]) {
                continue;
            }

            int prime_neighbors_including_self = 0;
            for (int dr = -1; dr <= 1; ++dr) {
                const int rr = rid + dr;
                const int rr_len = r_min + rr;
                for (int dc = -1; dc <= 1; ++dc) {
                    const int cc = col + dc;
                    if (cc < 0 || cc >= rr_len) {
                        continue;
                    }
                    prime_neighbors_including_self += row_prime[static_cast<std::size_t>(rr)][static_cast<std::size_t>(cc)] != 0;
                }
            }

            if (prime_neighbors_including_self >= 3) {
                for (int dr = -1; dr <= 1; ++dr) {
                    const int rr = rid + dr;
                    const int rr_len = r_min + rr;
                    for (int dc = -1; dc <= 1; ++dc) {
                        const int cc = col + dc;
                        if (cc < 0 || cc >= rr_len) {
                            continue;
                        }
                        if (row_prime[static_cast<std::size_t>(rr)][static_cast<std::size_t>(cc)]) {
                            row_good[static_cast<std::size_t>(rr)][static_cast<std::size_t>(cc)] = 1;
                        }
                    }
                }
            }
        }
    }

    u64 sum = 0;
    const int mid_idx = 2;
    const u64 mid_start = row_start(n);
    for (int col = 0; col < n; ++col) {
        if (row_good[static_cast<std::size_t>(mid_idx)][static_cast<std::size_t>(col)]) {
            sum += mid_start + static_cast<u64>(col);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve_row(8) != 60ULL) {
        std::cerr << "Checkpoint failed for S(8)" << '\n';
        return false;
    }
    if (solve_row(9) != 37ULL) {
        std::cerr << "Checkpoint failed for S(9)" << '\n';
        return false;
    }
    if (solve_row(10000) != 950007619ULL) {
        std::cerr << "Checkpoint failed for S(10000)" << '\n';
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

    if (options.single_row != 0) {
        std::cout << solve_row(options.single_row) << '\n';
    } else {
        std::cout << (solve_row(options.row_a) + solve_row(options.row_b)) << '\n';
    }
    return 0;
}
