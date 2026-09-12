#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using Wide = __uint128_t;

constexpr int MAX_BASE = 20;

Wide largest_value(const int a, const int b) {
    if (b >= 3 * a) return 0;
    if (b == 2 * a) return static_cast<Wide>(a) * (a - 1);

    std::vector<Wide> weights;
    Wide power_a = 1;
    Wide power_b = 1;
    Wide capacity = 0;

    // Let k be the first index with b^k >= 2*a^k; no higher digit is possible.
    while (power_b < 2 * power_a) {
        const Wide weight = 2 * power_a - power_b;
        capacity += (a - 1) * weight;
        weights.push_back(weight);
        power_a *= a;
        power_b *= b;
    }

    const Wide leading_weight = power_b - 2 * power_a;
    Wide value = std::min<Wide>(a - 1, capacity / leading_weight);
    Wide remaining = value * leading_weight;

    // w_i <= 1+(a-1)*sum_(j<i) w_j guarantees that descending greedy is exact.
    for (auto it = weights.rbegin(); it != weights.rend(); ++it) {
        const Wide digit = std::min<Wide>(a - 1, remaining / *it);
        remaining -= digit * *it;
        value = value * a + digit;
    }
    return value;
}

Wide sum_for_base(const int a) {
    Wide result = 0;
    for (int b = a + 1; b < 3 * a; ++b) result += largest_value(a, b);
    return result;
}

bool same_digits(Wide n, Wide doubled, const int a, const int b) {
    while (n != 0 || doubled != 0) {
        if (n % a != doubled % b) return false;
        n /= a;
        doubled /= b;
    }
    return true;
}

u64 brute_force(const int a, const int b) {
    u64 power_a = 1;
    u64 power_b = 1;
    while (power_b < 2 * power_a) {
        power_a *= a;
        power_b *= b;
    }
    u64 result = 0;
    for (u64 n = 1; n < power_a * a; ++n) {
        if (same_digits(n, 2 * n, a, b)) result = n;
    }
    return result;
}

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

void run_tests() {
    require(largest_value(3, 4) == 53, "F(3,4)");
    require(largest_value(9, 10) == 8'152'650, "F(9,10)");
    require(sum_for_base(3) == 72, "G(3)");
    for (int a = 2; a <= 7; ++a) {
        for (int b = a + 1; b <= 3 * a + 2; ++b) {
            require(largest_value(a, b) == brute_force(a, b),
                "exhaustive comparison for a=" + std::to_string(a) +
                ", b=" + std::to_string(b));
        }
    }
    for (int a = 2; a <= MAX_BASE; ++a) {
        for (int b = a + 1; b < 3 * a; ++b) {
            const Wide value = largest_value(a, b);
            require(same_digits(value, 2 * value, a, b), "matching base representations");
        }
    }
    std::cout << "All checks passed.\n";
}

std::string decimal(Wide value) {
    std::string result;
    do {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    } while (value != 0);
    std::reverse(result.begin(), result.end());
    return result;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
            return EXIT_SUCCESS;
        }
        if (argc != 1) throw std::invalid_argument("Usage: Euler1009 [--self-test]");
        Wide result = 0;
        for (int a = 2; a <= MAX_BASE; ++a) result += sum_for_base(a);
        std::cout << decimal(result) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
