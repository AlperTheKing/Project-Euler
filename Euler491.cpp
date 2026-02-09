#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

u128 factorial(const int n) {
    u128 v = 1;
    for (int i = 2; i <= n; ++i) {
        v *= static_cast<u128>(i);
    }
    return v;
}

u128 count_generic_pairs_div11(const int digits_count) {
    const int total_slots = 2 * digits_count;
    const int odd_slots = digits_count;
    const int even_slots = digits_count;

    std::vector<u128> fact(total_slots + 1, 1);
    for (int i = 2; i <= total_slots; ++i) {
        fact[static_cast<std::size_t>(i)] = fact[static_cast<std::size_t>(i - 1)] * static_cast<u128>(i);
    }

    std::vector<int> odd_take(static_cast<std::size_t>(digits_count), 0);
    u128 answer = 0;

    std::function<void(int, int, int)> dfs = [&](int digit, int used, int weighted) {
        if (digit == digits_count) {
            if (used != odd_slots) {
                return;
            }
            const int total_digit_sum = digits_count * (digits_count - 1) / 2;
            if ((weighted - total_digit_sum) % 11 != 0) {
                return;
            }

            u128 odd_all = fact[static_cast<std::size_t>(odd_slots)];
            u128 odd_zero_lead = 0;
            u128 even_all = fact[static_cast<std::size_t>(even_slots)];

            for (int d = 0; d < digits_count; ++d) {
                odd_all /= fact[static_cast<std::size_t>(odd_take[static_cast<std::size_t>(d)])];
                even_all /= fact[static_cast<std::size_t>(2 - odd_take[static_cast<std::size_t>(d)])];
            }
            if (odd_take[0] > 0) {
                odd_zero_lead = fact[static_cast<std::size_t>(odd_slots - 1)];
                odd_zero_lead /= fact[static_cast<std::size_t>(odd_take[0] - 1)];
                for (int d = 1; d < digits_count; ++d) {
                    odd_zero_lead /= fact[static_cast<std::size_t>(odd_take[static_cast<std::size_t>(d)])];
                }
            }

            const u128 odd_valid = odd_all - odd_zero_lead;
            answer += odd_valid * even_all;
            return;
        }

        for (int take = 0; take <= 2; ++take) {
            if (used + take > odd_slots) {
                continue;
            }
            odd_take[static_cast<std::size_t>(digit)] = take;
            dfs(digit + 1, used + take, weighted + digit * take);
        }
    };

    dfs(0, 0, 0);
    return answer;
}

u128 solve() {
    return count_generic_pairs_div11(10);
}

u128 brute_tiny(const int digits_count) {
    std::vector<int> digits;
    digits.reserve(static_cast<std::size_t>(2 * digits_count));
    for (int d = 0; d < digits_count; ++d) {
        digits.push_back(d);
        digits.push_back(d);
    }
    std::sort(digits.begin(), digits.end());

    u128 count = 0;
    do {
        if (digits[0] == 0) {
            continue;
        }
        int odd_sum = 0;
        int even_sum = 0;
        for (int i = 0; i < static_cast<int>(digits.size()); ++i) {
            if ((i % 2) == 0) {
                odd_sum += digits[static_cast<std::size_t>(i)];
            } else {
                even_sum += digits[static_cast<std::size_t>(i)];
            }
        }
        if ((odd_sum - even_sum) % 11 == 0) {
            ++count;
        }
    } while (std::next_permutation(digits.begin(), digits.end()));
    return count;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int d = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + d));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool run_checkpoints() {
    const u128 fast_small = count_generic_pairs_div11(4);
    const u128 brute_small = brute_tiny(4);
    if (fast_small != brute_small) {
        std::cerr << "Checkpoint failed: tiny brute-force cross-check (digits 0..3 twice)" << '\n';
        return false;
    }
    if (count_generic_pairs_div11(3) != brute_tiny(3)) {
        std::cerr << "Checkpoint failed: tiny brute-force cross-check (digits 0..2 twice)" << '\n';
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
    std::cout << to_string_u128(solve()) << '\n';
    return 0;
}
