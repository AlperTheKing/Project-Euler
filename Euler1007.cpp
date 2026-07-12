#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 MOD = 1'000'000'009LL;
constexpr int TARGET = 10'000'000;

i64 normalize(i64 value) {
    value %= MOD;
    return value < 0 ? value + MOD : value;
}

i64 multiply_mod(const i64 lhs, const i64 rhs) {
    return lhs * rhs % MOD;
}

i64 mod_pow(i64 base, i64 exponent) {
    i64 result = 1;
    while (exponent > 0) {
        if ((exponent & 1LL) != 0) {
            result = multiply_mod(result, base);
        }
        base = multiply_mod(base, base);
        exponent >>= 1LL;
    }
    return result;
}

struct RingElement {
    i64 constant;
    i64 root;
};

RingElement add(const RingElement lhs, const RingElement rhs) {
    return {normalize(lhs.constant + rhs.constant), normalize(lhs.root + rhs.root)};
}

RingElement subtract(const RingElement lhs, const RingElement rhs) {
    return {normalize(lhs.constant - rhs.constant), normalize(lhs.root - rhs.root)};
}

RingElement scale(const RingElement value, const i64 factor) {
    return {multiply_mod(value.constant, factor), multiply_mod(value.root, factor)};
}

// Represents a + b*r modulo r^2-r-1; the r coefficient supplies Fibonacci weights.
RingElement multiply_by_root(const RingElement value) {
    return {value.root, normalize(value.constant + value.root)};
}

RingElement multiply_by_root_plus_one(const RingElement value) {
    return {normalize(value.constant + value.root),
            normalize(value.constant + 2 * value.root)};
}

RingElement multiply_by_root_minus_one(const RingElement value) {
    return {normalize(value.root - value.constant), value.constant};
}

RingElement multiply_by_two_minus_root(const RingElement value) {
    return {normalize(2 * value.constant - value.root),
            normalize(value.root - value.constant)};
}

i64 alternating_sum(const int n) {
    const i64 inverse_four = mod_pow(4, MOD - 2);
    i64 a_scaled = 1;
    RingElement b_scaled{1, 0};
    RingElement c_older{0, 0};
    RingElement c_previous{1, 0};
    RingElement result_scaled{0, 0};
    i64 factorial = 1;

    for (i64 m = 1; m <= static_cast<i64>(n) + 1; ++m) {
        const i64 linear_factor = normalize(4 * m - 6);
        const i64 a_next = multiply_mod(linear_factor, a_scaled);
        const RingElement b_next = scale(multiply_by_root(b_scaled), linear_factor);

        RingElement c_next =
            scale(multiply_by_root_plus_one(c_previous), normalize(2 * (2 * m - 3)));
        const i64 quadratic_factor =
            multiply_mod(normalize(16 * normalize(m - 3)), normalize(m - 1));
        c_next = subtract(c_next, scale(multiply_by_root(c_older), quadratic_factor));

        RingElement numerator{normalize(-a_next), 0};
        numerator = add(numerator,
                        scale(multiply_by_root_minus_one({a_scaled, 0}), m % MOD));
        numerator = add(numerator, b_next);
        numerator = add(numerator,
                        scale(multiply_by_root_minus_one(b_scaled), m % MOD));
        numerator = subtract(numerator, c_next);
        if (m == 1) {
            numerator = add(numerator, {2, 2});
        }

        const RingElement adjusted =
            subtract(numerator,
                     scale(multiply_by_two_minus_root(result_scaled), (2 * m) % MOD));
        const RingElement result_next =
            scale(multiply_by_two_minus_root(adjusted), inverse_four);

        a_scaled = a_next;
        b_scaled = b_next;
        c_older = c_previous;
        c_previous = c_next;
        result_scaled = result_next;
        factorial = multiply_mod(factorial, m % MOD);
    }

    return multiply_mod(result_scaled.root, mod_pow(factorial, MOD - 2));
}

std::vector<i64> brute_values(const std::vector<i64>& fibonacci,
                              const int begin,
                              const int end) {
    if (end - begin == 1) {
        return {fibonacci[static_cast<std::size_t>(begin)]};
    }

    std::vector<i64> values;
    for (int split = begin + 1; split < end; ++split) {
        const std::vector<i64> left = brute_values(fibonacci, begin, split);
        const std::vector<i64> right = brute_values(fibonacci, split, end);
        for (const i64 lhs : left) {
            for (const i64 rhs : right) {
                values.push_back(lhs - rhs);
            }
        }
    }
    return values;
}

i64 brute_alternating_sum(const int n) {
    std::vector<i64> fibonacci(static_cast<std::size_t>(n + 1), 0);
    if (n >= 1) {
        fibonacci[1] = 1;
    }
    for (int index = 2; index <= n; ++index) {
        fibonacci[static_cast<std::size_t>(index)] =
            fibonacci[static_cast<std::size_t>(index - 1)] +
            fibonacci[static_cast<std::size_t>(index - 2)];
    }

    i64 total = 0;
    for (const i64 value : brute_values(fibonacci, 0, n + 1)) {
        total = normalize(total + value);
    }
    return total;
}

void require_checkpoint(const bool condition, const std::string& description) {
    if (!condition) {
        std::cerr << "Checkpoint failed: " << description << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void run_checkpoints() {
    for (int n = 0; n <= 10; ++n) {
        require_checkpoint(alternating_sum(n) == brute_alternating_sum(n),
                           "brute force comparison for n=" + std::to_string(n));
    }

    require_checkpoint(alternating_sum(3) == normalize(-6), "published A(3)");
    require_checkpoint(alternating_sum(10) == normalize(-177'666), "published A(10)");
    require_checkpoint(alternating_sum(100) == 71'792'794, "published A(100)");
}

}  // namespace

int main() {
    run_checkpoints();
    std::cout << alternating_sum(TARGET) << '\n';
    return 0;
}
