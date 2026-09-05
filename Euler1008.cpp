#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using Polynomial = std::array<u64, 10>;

constexpr u64 MOD = 1'000'000'007;
constexpr int TARGET = 10'000'000;
constexpr int BATCH_SIZE = 32'768;

u64 add(const u64 a, const u64 b) {
    const u64 sum = a + b;
    return sum >= MOD ? sum - MOD : sum;
}

u64 subtract(const u64 a, const u64 b) {
    return a >= b ? a - b : a + MOD - b;
}

u64 multiply(const u64 a, const u64 b) {
    return a * b % MOD;
}

u64 power(u64 base, u64 exponent) {
    u64 result = 1;
    while (exponent != 0) {
        if ((exponent & 1) != 0) result = multiply(result, base);
        base = multiply(base, base);
        exponent >>= 1;
    }
    return result;
}

struct Block {
    Polynomial product{1};
    Polynomial sum{};
    u64 factorial = 1;
};

template <bool FullSum>
Block calculate_block(const int begin, const int end) {
    Block result;
    std::vector<u64> inverse(BATCH_SIZE + 1);

    for (int first = begin; first < end; first += BATCH_SIZE) {
        const int count = std::min(BATCH_SIZE, end - first);
        inverse[0] = 1;
        for (int i = 1; i <= count; ++i) {
            const u64 k = static_cast<u64>(first + i - 1);
            inverse[i] = multiply(inverse[i - 1], multiply(k, 2 * k - 1));
        }
        u64 suffix_inverse = power(inverse[count], MOD - 2);
        for (int i = count; i >= 1; --i) {
            const u64 k = static_cast<u64>(first + i - 1);
            const u64 weight = multiply(suffix_inverse, inverse[i - 1]);
            suffix_inverse = multiply(suffix_inverse, multiply(k, 2 * k - 1));
            inverse[i] = weight;
        }

        for (int i = 1; i <= count; ++i) {
            const u64 k = static_cast<u64>(first + i - 1);
            const u64 weight = inverse[i];
            const u64 inverse_k = multiply(2 * k - 1, weight);
            const u64 inverse_square = multiply(inverse_k, inverse_k);
            result.factorial = multiply(result.factorial, k);

            // P_N/x = sum A_(k-1)/(k(2k-1)), A_k = product_(j<=k)(1-x/j^2).
            if constexpr (FullSum) {
                for (int degree = 9; degree >= 1; --degree) {
                    result.sum[degree] = add(result.sum[degree],
                                             multiply(weight, result.product[degree]));
                    result.product[degree] = subtract(result.product[degree],
                        multiply(inverse_square, result.product[degree - 1]));
                }
                result.sum[0] = add(result.sum[0], weight);
            } else {
                result.sum[9] = add(result.sum[9], multiply(weight, result.product[9]));
                for (int degree = 9; degree >= 1; --degree) {
                    result.product[degree] = subtract(result.product[degree],
                        multiply(inverse_square, result.product[degree - 1]));
                }
            }
        }
    }
    return result;
}

Polynomial product(const Polynomial& a, const Polynomial& b) {
    Polynomial result{};
    for (int i = 0; i <= 9; ++i) {
        for (int j = 0; j <= i; ++j) {
            result[i] = add(result[i], multiply(a[j], b[i - j]));
        }
    }
    return result;
}

u64 solve(const int n, unsigned thread_count) {
    thread_count = std::max(1U, std::min(thread_count, static_cast<unsigned>(n)));
    Block total;
    if (thread_count == 1) {
        total = calculate_block<false>(1, n + 1);
    } else {
        std::vector<Block> blocks(thread_count);
        std::vector<std::thread> workers;
        for (unsigned index = 0; index < thread_count; ++index) {
            const int begin = 1 + static_cast<int>(static_cast<u64>(n) * index / thread_count);
            const int end = 1 + static_cast<int>(static_cast<u64>(n) * (index + 1) / thread_count);
            workers.emplace_back([&, index, begin, end] {
                blocks[index] = calculate_block<true>(begin, end);
            });
        }
        for (auto& worker : workers) worker.join();
        for (const Block& block : blocks) {
            const Polynomial contribution = product(total.product, block.sum);
            total.sum[9] = add(total.sum[9], contribution[9]);
            total.product = product(total.product, block.product);
            total.factorial = multiply(total.factorial, block.factorial);
        }
    }

    const u64 factorial_square = multiply(total.factorial, total.factorial);
    u64 leading = multiply(static_cast<u64>(n),
        power(multiply(2ULL * n - 1, factorial_square), MOD - 2));
    if (n % 2 == 0) leading = subtract(0, leading);
    if (leading == 1) return total.sum[9];

    // Otherwise the minimum monic polynomial is P_N + x*product_(k<=N)(x-k^2).
    const u64 scale = n % 2 == 0 ? factorial_square : subtract(0, factorial_square);
    return add(total.sum[9], multiply(scale, total.product[9]));
}

std::vector<u64> append_root(const std::vector<u64>& polynomial, const u64 root) {
    std::vector<u64> result(polynomial.size() + 1);
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        result[i] = subtract(result[i], multiply(root, polynomial[i]));
        result[i + 1] = add(result[i + 1], polynomial[i]);
    }
    return result;
}

u64 interpolate_directly(const int n) {
    std::vector<u64> polynomial(static_cast<std::size_t>(n + 2));
    for (int i = 1; i <= n; ++i) {
        std::vector<u64> basis{1};
        u64 denominator = 1;
        for (int j = 0; j <= n; ++j) {
            if (j == i) continue;
            basis = append_root(basis, static_cast<u64>(j) * j);
            denominator = multiply(denominator,
                subtract(static_cast<u64>(i) * i, static_cast<u64>(j) * j));
        }
        const u64 scale = multiply(static_cast<u64>(i), power(denominator, MOD - 2));
        for (std::size_t degree = 0; degree < basis.size(); ++degree) {
            polynomial[degree] = add(polynomial[degree], multiply(scale, basis[degree]));
        }
    }
    if (polynomial[n] != 1) {
        std::vector<u64> vanishing{1};
        for (int i = 0; i <= n; ++i) {
            vanishing = append_root(vanishing, static_cast<u64>(i) * i);
        }
        for (std::size_t degree = 0; degree < vanishing.size(); ++degree) {
            polynomial[degree] = add(polynomial[degree], vanishing[degree]);
        }
    }
    return polynomial.size() > 10 ? polynomial[10] : 0;
}

void run_tests() {
    for (int n = 1; n <= 40; ++n) {
        const u64 expected = interpolate_directly(n);
        if (solve(n, 1) != expected || solve(n, 3) != expected) {
            throw std::runtime_error("Interpolation check failed for n=" + std::to_string(n));
        }
    }
    for (const int n : {BATCH_SIZE - 1, BATCH_SIZE, BATCH_SIZE + 1, 2 * BATCH_SIZE + 3}) {
        if (solve(n, 1) != solve(n, 3)) {
            throw std::runtime_error("Batch boundary check failed for n=" + std::to_string(n));
        }
    }
    std::cout << "All checks passed.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        unsigned thread_count = std::min(16U, std::max(1U, std::thread::hardware_concurrency()));
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "--threads") {
            const std::string value = argv[2];
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos) {
                throw std::invalid_argument("Thread count must be an integer from 1 to 128.");
            }
            const unsigned long parsed = std::stoul(value);
            if (parsed < 1 || parsed > 128) {
                throw std::invalid_argument("Thread count must be an integer from 1 to 128.");
            }
            thread_count = static_cast<unsigned>(parsed);
        } else if (argc != 1) {
            throw std::invalid_argument("Usage: Euler1008 [--threads COUNT | --self-test]");
        }
        std::cout << solve(TARGET, thread_count) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
