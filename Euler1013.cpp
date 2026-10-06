#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
constexpr u64 LIMIT = 1'000'000'000'000'000;

constexpr u64 fourth_power(const u64 n) {
    const u64 square = n * n;
    return square * square;
}

unsigned fourth_root(const u64 n) {
    unsigned low = 0, high = 65'536;
    while (high - low > 1) {
        const unsigned middle = low + (high - low) / 2;
        if (fourth_power(middle) <= n) low = middle;
        else high = middle;
    }
    return low;
}

u64 solve(const u64 limit) {
    if (limit > LIMIT) throw std::invalid_argument("Supported limit is at most 10^15");
    if (limit <= 1) return 0;
    const unsigned bound = fourth_root(limit - 1);
    std::vector<bool> composite(bound + 1);
    std::vector<u64> primes;
    for (unsigned p = 2; p <= bound; ++p) {
        if (composite[p]) continue;
        if (p % 5 == 1) primes.push_back(p);
        for (unsigned multiple = p * p; multiple <= bound; multiple += p) composite[multiple] = true;
    }
    u64 sum = 0;
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const u64 p = primes[i];
        if (p > bound / p) break;
        for (std::size_t j = i + 1; j < primes.size(); ++j) {
            const u64 q = primes[j];
            if (q > bound / p) break;
            sum += fourth_power(p * q);
        }
    }
    return sum;
}

bool qualifies(u64 n) {
    u64 count = 1, sum = 1;
    for (u64 p = 2; p <= n / p; ++p) {
        if (n % p != 0) continue;
        unsigned exponent = 0;
        u64 power = 1, factor_sum = 1;
        do {
            n /= p;
            ++exponent;
            power *= p;
            factor_sum += power;
        } while (n % p == 0);
        count *= exponent + 1;
        sum *= factor_sum;
    }
    if (n > 1) {
        count *= 2;
        sum *= n + 1;
    }
    return count == 25 && sum % 25 == 0;
}

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

void run_tests() {
    constexpr u64 FIRST = 13'521'270'961;
    require(qualifies(FIRST), "given first number");
    require(solve(FIRST) == 0 && solve(FIRST + 1) == FIRST, "first number and strict cutoff");
    require(solve(0) == 0 && solve(1) == 0, "empty ranges");
    u64 expected = 0;
    unsigned count = 0;
    for (u64 base = 1; fourth_power(base) < LIMIT; ++base) {
        const u64 n = fourth_power(base);
        if (!qualifies(n)) continue;
        require(solve(n) == expected, "exclusive cutoff at each qualifying number");
        expected += n;
        ++count;
        require(solve(n + 1) == expected, "inclusive successor at each qualifying number");
    }
    require(solve(LIMIT) == expected, "full-range factorization check");
    std::cout << "All checks passed; " << count << " qualifying numbers.\n";
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
            return EXIT_SUCCESS;
        }
        if (argc != 1) throw std::invalid_argument("Usage: Euler1013 [--self-test]");
        std::cout << solve(LIMIT) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
