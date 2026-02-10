#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

// Project Euler 571: Super Pandigital Numbers
//
// An n-super-pandigital number must be pandigital in base n, so its base-n
// representation contains all digits 0..n-1 at least once. Hence it has at
// least n digits in base n, and the smallest candidates are exactly the
// n-digit base-n pandigitals, i.e. permutations of 0..n-1 with a nonzero
// leading digit (each digit appears exactly once).
//
// We enumerate those permutations in increasing numeric order (lexicographic
// order of the base-n digit sequence) and test pandigitality in bases 2..n-1.
// The first k hits are the k smallest n-super-pandigital numbers.

using u64 = std::uint64_t;

static inline bool is_pandigital_in_base(u64 x, int base) {
    const std::uint32_t all = (base == 32) ? 0xFFFFFFFFu : ((1u << base) - 1u);
    std::uint32_t mask = 0;
    while (x) {
        mask |= 1u << (x % static_cast<u64>(base));
        if (mask == all) return true;  // early success
        x /= static_cast<u64>(base);
    }
    return mask == all;
}

static inline bool is_super_pandigital(u64 x, int n) {
    // Base n is guaranteed by construction in the generator below.
    for (int b = n - 1; b >= 2; --b) {
        if (!is_pandigital_in_base(x, b)) return false;
    }
    return true;
}

static std::vector<u64> k_smallest_super_pandigital(int n, int k) {
    std::vector<int> digits(n);
    std::iota(digits.begin(), digits.end(), 0);
    std::swap(digits[0], digits[1]);  // 1 0 2 3 ... is the smallest with nonzero leading digit.

    std::vector<u64> out;
    out.reserve(k);

    do {
        // Convert base-n digits -> value.
        u64 v = 0;
        for (int d : digits) v = v * static_cast<u64>(n) + static_cast<u64>(d);

        if (is_super_pandigital(v, n)) {
            out.push_back(v);
            if (static_cast<int>(out.size()) == k) break;
        }
    } while (std::next_permutation(digits.begin(), digits.end()));

    return out;
}

int main() {
    // Validation from the statement.
    {
        const auto v5 = k_smallest_super_pandigital(5, 1);
        if (v5.size() != 1 || v5[0] != 978ULL) {
            std::cerr << "Validation failed for smallest 5-super-pandigital\n";
            return 1;
        }
        const auto v10 = k_smallest_super_pandigital(10, 10);
        if (v10.size() != 10 || v10[0] != 1093265784ULL) {
            std::cerr << "Validation failed for smallest 10-super-pandigital\n";
            return 1;
        }
        const u64 sum10 = std::accumulate(v10.begin(), v10.end(), 0ULL);
        if (sum10 != 20319792309ULL) {
            std::cerr << "Validation failed for sum of 10 smallest 10-super-pandigital\n";
            return 1;
        }
    }

    const auto v12 = k_smallest_super_pandigital(12, 10);
    const u64 ans = std::accumulate(v12.begin(), v12.end(), 0ULL);
    std::cout << ans << "\n";
    return 0;
}

