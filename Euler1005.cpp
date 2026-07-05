#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;
using i64 = std::int64_t;

constexpr int kDefaultTarget = 2026;
constexpr i64 kMod = 1'000'000'000LL;

struct Options {
    int target = kDefaultTarget;
    bool run_checkpoints = true;
};

struct PrimeDp {
    int target = 0;
    std::vector<int> primes;
    std::vector<std::vector<cpp_int>> suffix_count;
};

std::vector<int> primes_up_to(const int limit) {
    if (limit < 2) {
        return {};
    }

    std::vector<bool> composite(static_cast<std::size_t>(limit + 1), false);
    for (int p = 2; p * p <= limit; ++p) {
        if (!composite[static_cast<std::size_t>(p)]) {
            for (int multiple = p * p; multiple <= limit; multiple += p) {
                composite[static_cast<std::size_t>(multiple)] = true;
            }
        }
    }

    std::vector<int> primes;
    for (int value = 2; value <= limit; ++value) {
        if (!composite[static_cast<std::size_t>(value)]) {
            primes.push_back(value);
        }
    }
    return primes;
}

PrimeDp build_dp(const int target) {
    PrimeDp dp;
    dp.target = target;
    dp.primes = primes_up_to(target);
    const int count = static_cast<int>(dp.primes.size());
    dp.suffix_count.assign(static_cast<std::size_t>(count + 1),
                           std::vector<cpp_int>(static_cast<std::size_t>(target + 1), 0));
    dp.suffix_count[static_cast<std::size_t>(count)][0] = 1;

    for (int i = count - 1; i >= 0; --i) {
        const int prime = dp.primes[static_cast<std::size_t>(i)];
        for (int sum = 0; sum <= target; ++sum) {
            cpp_int ways = dp.suffix_count[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(sum)];
            if (sum >= prime) {
                ways += dp.suffix_count[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(sum - prime)];
            }
            dp.suffix_count[static_cast<std::size_t>(i)][static_cast<std::size_t>(sum)] = ways;
        }
    }

    return dp;
}

std::vector<int> kth_prime_list(const PrimeDp& dp, cpp_int rank) {
    if (rank <= 0 || rank > dp.suffix_count[0][static_cast<std::size_t>(dp.target)]) {
        throw std::runtime_error("Requested rank is outside the available prime lists");
    }

    int remaining = dp.target;
    int start = 0;
    std::vector<int> result;

    while (remaining > 0) {
        bool found = false;
        for (int i = start; i < static_cast<int>(dp.primes.size()); ++i) {
            const int prime = dp.primes[static_cast<std::size_t>(i)];
            if (prime > remaining) {
                break;
            }

            const cpp_int& block =
                dp.suffix_count[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(remaining - prime)];
            if (rank > block) {
                rank -= block;
                continue;
            }

            result.push_back(prime);
            remaining -= prime;
            start = i + 1;
            found = true;
            break;
        }

        if (!found) {
            throw std::runtime_error("Could not unrank the requested prime list");
        }
    }

    return result;
}

std::vector<int> median_prime_list(const PrimeDp& dp) {
    const cpp_int total = dp.suffix_count[0][static_cast<std::size_t>(dp.target)];
    if (total == 0) {
        throw std::runtime_error("There is no prime list for the requested target");
    }
    return kth_prime_list(dp, (total + 1) / 2);
}

i64 product_mod(const std::vector<int>& values) {
    i64 product = 1;
    for (const int value : values) {
        product = product * value % kMod;
    }
    return product;
}

void enumerate_bruteforce(const std::vector<int>& primes,
                          const int start,
                          const int remaining,
                          std::vector<int>& current,
                          std::vector<std::vector<int>>& lists) {
    if (remaining == 0) {
        lists.push_back(current);
        return;
    }

    for (int i = start; i < static_cast<int>(primes.size()); ++i) {
        const int prime = primes[static_cast<std::size_t>(i)];
        if (prime > remaining) {
            break;
        }
        current.push_back(prime);
        enumerate_bruteforce(primes, i + 1, remaining - prime, current, lists);
        current.pop_back();
    }
}

std::vector<std::vector<int>> brute_lists(const int target) {
    std::vector<std::vector<int>> lists;
    std::vector<int> current;
    const std::vector<int> primes = primes_up_to(target);
    enumerate_bruteforce(primes, 0, target, current, lists);
    return lists;
}

void require_checkpoint(const bool ok, const std::string& message) {
    if (!ok) {
        std::cerr << "Checkpoint failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void run_checkpoints() {
    {
        const PrimeDp dp = build_dp(20);
        const std::vector<std::vector<int>> expected{{2, 5, 13}, {2, 7, 11}, {3, 17}, {7, 13}};
        require_checkpoint(brute_lists(20) == expected, "lexicographic lists for 20");
        require_checkpoint(dp.suffix_count[0][20] == 4, "count for 20");
        require_checkpoint(median_prime_list(dp) == std::vector<int>({2, 7, 11}), "median list for 20");
        require_checkpoint(product_mod(median_prime_list(dp)) == 154, "product for 20");
    }

    for (int target = 2; target <= 60; ++target) {
        const PrimeDp dp = build_dp(target);
        const std::vector<std::vector<int>> lists = brute_lists(target);
        require_checkpoint(dp.suffix_count[0][static_cast<std::size_t>(target)] == lists.size(),
                           "DP count matches brute force for " + std::to_string(target));

        if (!lists.empty()) {
            const std::vector<int> median = median_prime_list(dp);
            require_checkpoint(median == lists[(lists.size() - 1) / 2],
                               "median matches brute force for " + std::to_string(target));

            const std::vector<int> first = kth_prime_list(dp, 1);
            const std::vector<int> middle = kth_prime_list(dp, cpp_int((lists.size() + 1) / 2));
            const std::vector<int> last = kth_prime_list(dp, cpp_int(lists.size()));
            require_checkpoint(first == lists.front(), "first rank for " + std::to_string(target));
            require_checkpoint(middle == lists[(lists.size() - 1) / 2],
                               "middle rank for " + std::to_string(target));
            require_checkpoint(last == lists.back(), "last rank for " + std::to_string(target));
        }
    }

    std::cerr << "Validation checkpoints passed.\n";
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << '\n';
            return false;
        }
        options.target = std::stoi(arg);
    }

    if (options.target < 0) {
        std::cerr << "Target must be nonnegative.\n";
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

    if (options.run_checkpoints) {
        run_checkpoints();
    }

    try {
        const PrimeDp dp = build_dp(options.target);
        const std::vector<int> median = median_prime_list(dp);
        std::cout << std::setw(9) << std::setfill('0') << product_mod(median) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 2;
    }

    return 0;
}
