#include <array>
#include <bit>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 MOD = 1'000'000'007;
constexpr u64 TARGET = 1'000'000'000'000'000'000;
constexpr u64 MAX_LIMIT = std::numeric_limits<u64>::max() - 1;

void require(const bool condition, const char* description) {
    if (!condition) throw std::runtime_error(std::string("Check failed: ") + description);
}

struct Modular {
    using Value = u64;

    static Value add(const Value a, const Value b) {
        const Value sum = a + b;
        return sum >= MOD ? sum - MOD : sum;
    }

    static Value multiply(const Value a, const Value b) {
        return static_cast<Value>(static_cast<u128>(a) * b % MOD);
    }
};

struct Exact {
    using Value = u128;

    static Value add(const Value a, const Value b) {
        Value sum;
        require(!__builtin_add_overflow(a, b, &sum), "exact sum fits in 128 bits");
        return sum;
    }

    static Value multiply(const Value a, const Value b) {
        Value product;
        require(!__builtin_mul_overflow(a, b, &product), "exact product fits in 128 bits");
        return product;
    }
};

template <typename Ring>
struct Matrix {
    using Value = typename Ring::Value;
    std::array<std::array<Value, 2>, 2> entry;

    static Matrix zero() { return {{{{0, 0}, {0, 0}}}}; }
    static Matrix identity() { return {{{{1, 0}, {0, 1}}}}; }
    static Matrix fibonacci() { return {{{{1, 1}, {1, 0}}}}; }

    Matrix operator+(const Matrix& other) const {
        Matrix result{};
        for (unsigned i = 0; i < 2; ++i) {
            for (unsigned j = 0; j < 2; ++j) result.entry[i][j] = Ring::add(entry[i][j], other.entry[i][j]);
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result{};
        for (unsigned i = 0; i < 2; ++i) {
            for (unsigned j = 0; j < 2; ++j) {
                result.entry[i][j] = Ring::add(Ring::multiply(entry[i][0], other.entry[0][j]),
                                               Ring::multiply(entry[i][1], other.entry[1][j]));
            }
        }
        return result;
    }
};

template <typename Ring>
Matrix<Ring> power(Matrix<Ring> base, u64 exponent) {
    auto result = Matrix<Ring>::identity();
    while (exponent != 0) {
        if ((exponent & 1) != 0) result = result * base;
        exponent >>= 1;
        if (exponent != 0) base = base * base;
    }
    return result;
}

template <typename Ring>
typename Ring::Value modified_sum(const u64 limit) {
    using M = Matrix<Ring>;
    require(limit <= MAX_LIMIT, "N + 1 fits in 64 bits");
    const u64 count = limit + 1;
    const unsigned bits = static_cast<unsigned>(std::bit_width(count));
    std::vector<M> factor, block{M::identity()};
    M square = M::fibonacci();
    for (unsigned i = 0; i < bits; ++i) {
        const M shifted = square * M::fibonacci();
        block.push_back(block.back() * (M::identity() + shifted));
        factor.push_back(shifted);
        if (i + 1 < bits) square = square * square;
    }
    M total = M::zero();
    M prefix = M::identity();
    for (unsigned bit = bits; bit-- > 0;) {
        if (((count >> bit) & 1) != 0) {
            total = total + prefix * block[bit];
            prefix = prefix * factor[bit];
        }
    }
    return total.entry[0][1];
}

struct Quadratic {
    u64 a, b;
};

Quadratic operator+(const Quadratic& x, const Quadratic& y) {
    return {Modular::add(x.a, y.a), Modular::add(x.b, y.b)};
}

Quadratic operator*(const Quadratic& x, const Quadratic& y) {
    const u128 a = static_cast<u128>(x.a) * y.a + 5 * (static_cast<u128>(x.b) * y.b);
    const u128 b = static_cast<u128>(x.a) * y.b + static_cast<u128>(x.b) * y.a;
    return {static_cast<u64>(a % MOD), static_cast<u64>(b % MOD)};
}

Quadratic power(Quadratic base, u64 exponent) {
    Quadratic result{1, 0};
    while (exponent != 0) {
        if ((exponent & 1) != 0) result = result * base;
        exponent >>= 1;
        if (exponent != 0) base = base * base;
    }
    return result;
}

Quadratic weighted_sum(const u64 limit, const Quadratic& base, const Quadratic& phi) {
    if (limit == 0) return {1, 0};
    if ((limit & 1) != 0) return (Quadratic{1, 0} + base * phi) * weighted_sum(limit / 2, base * base, phi);
    return weighted_sum(limit - 1, base, phi) + power(base, limit) * power(phi, std::popcount(limit));
}

u64 binet_sum(const u64 limit) {
    constexpr u64 HALF = (MOD + 1) / 2;
    const Quadratic phi{HALF, HALF};
    const Quadratic total = weighted_sum(limit, phi, phi);
    return Modular::add(total.b, total.b);
}

template <typename Ring>
std::vector<typename Ring::Value> brute_force_prefixes(const unsigned limit) {
    std::vector<typename Ring::Value> fibonacci(limit + 65), prefixes;
    fibonacci[1] = 1;
    for (std::size_t i = 2; i < fibonacci.size(); ++i) fibonacci[i] = Ring::add(fibonacci[i - 1], fibonacci[i - 2]);
    typename Ring::Value total = 0;
    for (unsigned k = 0; k <= limit; ++k) {
        total = Ring::add(total, fibonacci[k + static_cast<unsigned>(std::popcount(k))]);
        prefixes.push_back(total);
    }
    return prefixes;
}

void run_tests() {
    require(modified_sum<Exact>(10) == 369, "given S(10), matrix product");
    require(modified_sum<Modular>(10) == 369, "given S(10) modulo the prime");
    require(binet_sum(10) == 369, "given S(10), Binet method");
    const u128 s50 = 161'052'232'733;
    require(modified_sum<Exact>(50) == s50, "given S(50), matrix product");
    require(modified_sum<Modular>(50) == s50 % MOD, "given S(50) modulo the prime");
    require(binet_sum(50) == s50 % MOD, "given S(50), Binet method");

    const auto exact = brute_force_prefixes<Exact>(100);
    require(exact[10] == 369 && exact[50] == s50, "given values by direct summation");
    for (unsigned limit = 0; limit <= 100; ++limit) {
        require(modified_sum<Exact>(limit) == exact[limit], "exact sums for every N up to 100");
    }
    const auto prefixes = brute_force_prefixes<Modular>(1'000'000);
    for (unsigned limit = 0; limit <= 4096; ++limit) {
        require(modified_sum<Modular>(limit) == prefixes[limit], "matrix product for every N up to 4096");
        require(binet_sum(limit) == prefixes[limit], "Binet method for every N up to 4096");
    }
    for (const unsigned limit : {0u, 1u, 2u, 7u, 8u, 15u, 16u, 31u, 32u, 1000u, 123'456u, 1'000'000u}) {
        require(modified_sum<Modular>(limit) == prefixes[limit], "matrix product against direct summation");
        require(binet_sum(limit) == prefixes[limit], "Binet method against direct summation");
    }

    const u64 answer = modified_sum<Modular>(TARGET);
    require(answer == binet_sum(TARGET), "both methods at N = 10^18");
    const u64 step = Modular::add(answer, MOD - modified_sum<Modular>(TARGET - 1));
    const u64 last = power(Matrix<Modular>::fibonacci(), TARGET + std::popcount(TARGET)).entry[0][1];
    require(step == last, "S(N) - S(N - 1) = F(N + s(N)) at N = 10^18");
    require(modified_sum<Modular>(MAX_LIMIT) == binet_sum(MAX_LIMIT), "both methods at the largest N");
    std::cout << "All checks passed: S(10) and S(50), exact sums for N up to 100, direct sums for every N up to 4096 "
                 "and at 12 values up to 1000000, the Binet method, the last term at N = 10^18 and both methods at the "
                 "largest N.\n";
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
        } else {
            if (argc > 2) throw std::invalid_argument("Usage: Euler1020 [N | --self-test]");
            u64 limit = TARGET;
            if (argc == 2) {
                const std::string argument(argv[1]);
                const auto parsed = std::from_chars(argument.data(), argument.data() + argument.size(), limit);
                if (parsed.ec != std::errc{} || parsed.ptr != argument.data() + argument.size() || limit > MAX_LIMIT) {
                    throw std::invalid_argument("N must be an integer from 0 to 18446744073709551614");
                }
            }
            std::cout << modified_sum<Modular>(limit) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
