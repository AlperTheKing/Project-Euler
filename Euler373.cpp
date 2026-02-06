#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int kTargetRadius = 10'000'000;

struct U192 {
    u64 lo = 0;
    u64 mid = 0;
    u64 hi = 0;

    static U192 one() {
        U192 value;
        value.lo = 1;
        return value;
    }

    void multiply_small(u64 factor) {
        const u128 t0 = static_cast<u128>(lo) * factor;
        lo = static_cast<u64>(t0);
        u64 carry = static_cast<u64>(t0 >> 64);

        const u128 t1 = static_cast<u128>(mid) * factor + carry;
        mid = static_cast<u64>(t1);
        carry = static_cast<u64>(t1 >> 64);

        const u128 t2 = static_cast<u128>(hi) * factor + carry;
        hi = static_cast<u64>(t2);
    }

    bool operator==(const U192& other) const {
        return lo == other.lo && mid == other.mid && hi == other.hi;
    }
};

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value > 0) {
        digits.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

bool circumradius_equals(int a, int b, int c, int radius) {
    // Exact check of:
    // a^2 b^2 c^2 == R^2 (a+b+c)(a+b-c)(a-b+c)(-a+b+c)
    U192 lhs = U192::one();
    lhs.multiply_small(static_cast<u64>(radius));
    lhs.multiply_small(static_cast<u64>(radius));
    lhs.multiply_small(static_cast<u64>(a + b + c));
    lhs.multiply_small(static_cast<u64>(a + b - c));
    lhs.multiply_small(static_cast<u64>(a - b + c));
    lhs.multiply_small(static_cast<u64>(-a + b + c));

    U192 rhs = U192::one();
    rhs.multiply_small(static_cast<u64>(a));
    rhs.multiply_small(static_cast<u64>(a));
    rhs.multiply_small(static_cast<u64>(b));
    rhs.multiply_small(static_cast<u64>(b));
    rhs.multiply_small(static_cast<u64>(c));
    rhs.multiply_small(static_cast<u64>(c));

    return lhs == rhs;
}

class Euler373Solver {
public:
    explicit Euler373Solver(int radius_limit)
        : radius_limit_(radius_limit),
          smallest_prime_factor_(static_cast<std::size_t>(radius_limit) + 1),
          denominator_to_numerators_(static_cast<std::size_t>(radius_limit) + 1) {
        build_smallest_prime_factor();
        build_denominator_table();
    }

    u128 solve(bool allow_multithreading, unsigned requested_threads = 0) const {
        return sum_up_to(radius_limit_, allow_multithreading, requested_threads);
    }

    u128 sum_up_to(int max_radius,
                   bool allow_multithreading,
                   unsigned requested_threads = 0) const {
        if (max_radius > radius_limit_) {
            return 0;
        }

        unsigned threads = requested_threads;
        if (threads == 0) {
            threads = std::thread::hardware_concurrency();
            if (threads == 0) {
                threads = 1;
            }
        }

        if (!allow_multithreading || threads <= 1 || max_radius < 50'000) {
            threads = 1;
        } else {
            threads = std::min<unsigned>(threads, static_cast<unsigned>(max_radius));
        }

        if (threads == 1) {
            std::vector<int> divisors;
            std::vector<int> sides;
            u128 total = 0;
            for (int radius = 1; radius <= max_radius; ++radius) {
                total += static_cast<u128>(radius) *
                         count_triangles_for_radius(radius, divisors, sides);
            }
            return total;
        }

        std::atomic<int> next_radius{1};
        std::vector<u128> partial(threads, 0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                std::vector<int> divisors;
                std::vector<int> sides;
                u128 local = 0;

                while (true) {
                    const int radius = next_radius.fetch_add(1, std::memory_order_relaxed);
                    if (radius > max_radius) {
                        break;
                    }

                    local += static_cast<u128>(radius) *
                             count_triangles_for_radius(radius, divisors, sides);
                }

                partial[tid] = local;
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        u128 total = 0;
        for (const u128 value : partial) {
            total += value;
        }
        return total;
    }

    u64 count_triangles_for_radius(int radius) const {
        std::vector<int> divisors;
        std::vector<int> sides;
        return count_triangles_for_radius(radius, divisors, sides);
    }

private:
    int radius_limit_ = 0;
    std::vector<int> smallest_prime_factor_;
    std::vector<std::vector<int>> denominator_to_numerators_;

    void build_smallest_prime_factor() {
        for (int value = 0; value <= radius_limit_; ++value) {
            smallest_prime_factor_[static_cast<std::size_t>(value)] = value;
        }
        if (radius_limit_ >= 1) {
            smallest_prime_factor_[1] = 1;
        }

        for (int p = 2; static_cast<std::int64_t>(p) * p <= radius_limit_; ++p) {
            if (smallest_prime_factor_[static_cast<std::size_t>(p)] != p) {
                continue;
            }

            for (int multiple = p * p; multiple <= radius_limit_; multiple += p) {
                if (smallest_prime_factor_[static_cast<std::size_t>(multiple)] == multiple) {
                    smallest_prime_factor_[static_cast<std::size_t>(multiple)] = p;
                }
            }
        }
    }

    void build_denominator_table() {
        const int m_limit = static_cast<int>(std::sqrt(static_cast<long double>(radius_limit_)));

        for (int m = 2; m <= m_limit; ++m) {
            for (int n = 1; n < m; ++n) {
                if (((m + n) & 1) == 0) {
                    continue;
                }
                if (std::gcd(m, n) != 1) {
                    continue;
                }

                const int denom = m * m + n * n;
                if (denom > radius_limit_) {
                    break;
                }

                auto& numerators = denominator_to_numerators_[static_cast<std::size_t>(denom)];
                numerators.push_back(m * m - n * n);
                numerators.push_back(2 * m * n);
            }
        }

        for (int denom = 1; denom <= radius_limit_; ++denom) {
            auto& numerators = denominator_to_numerators_[static_cast<std::size_t>(denom)];
            if (numerators.empty()) {
                continue;
            }
            std::sort(numerators.begin(), numerators.end());
            numerators.erase(std::unique(numerators.begin(), numerators.end()),
                             numerators.end());
        }
    }

    void divisors_of(int value, std::vector<int>& out) const {
        out.clear();
        out.push_back(1);

        while (value > 1) {
            const int prime = smallest_prime_factor_[static_cast<std::size_t>(value)];
            int exponent = 0;
            while (value % prime == 0) {
                value /= prime;
                ++exponent;
            }

            const std::size_t base_size = out.size();
            int multiplier = 1;
            for (int i = 1; i <= exponent; ++i) {
                multiplier *= prime;
                for (std::size_t idx = 0; idx < base_size; ++idx) {
                    out.push_back(out[idx] * multiplier);
                }
            }
        }
    }

    u64 count_triangles_for_radius(int radius,
                                   std::vector<int>& divisors,
                                   std::vector<int>& sides) const {
        divisors_of(radius, divisors);

        sides.clear();
        sides.push_back(2 * radius); // Diameter (sin(theta)=1).

        for (const int denom : divisors) {
            const auto& numerators =
                denominator_to_numerators_[static_cast<std::size_t>(denom)];
            if (numerators.empty()) {
                continue;
            }

            const int scale = 2 * (radius / denom);
            for (const int numer : numerators) {
                sides.push_back(scale * numer);
            }
        }

        std::sort(sides.begin(), sides.end());
        sides.erase(std::unique(sides.begin(), sides.end()), sides.end());

        const int count = static_cast<int>(sides.size());
        u64 triangles = 0;

        for (int i = 0; i < count; ++i) {
            const int a = sides[static_cast<std::size_t>(i)];
            for (int j = i; j < count; ++j) {
                const int b = sides[static_cast<std::size_t>(j)];
                const int ab = a + b;
                for (int k = j; k < count; ++k) {
                    const int c = sides[static_cast<std::size_t>(k)];
                    if (ab <= c) {
                        break;
                    }

                    if (circumradius_equals(a, b, c, radius)) {
                        ++triangles;
                    }
                }
            }
        }

        return triangles;
    }
};

u64 brute_force_triangle_count_for_radius(int radius) {
    const int max_side = 2 * radius;
    u64 count = 0;

    for (int a = 1; a <= max_side; ++a) {
        for (int b = a; b <= max_side; ++b) {
            const int c_max = std::min(max_side, a + b - 1);
            for (int c = b; c <= c_max; ++c) {
                if (circumradius_equals(a, b, c, radius)) {
                    ++count;
                }
            }
        }
    }

    return count;
}

bool run_validation_checkpoints() {
    {
        // Structural validation: compare the optimized candidate enumeration
        // against brute force on small radii.
        Euler373Solver solver(60);
        for (int radius = 1; radius <= 60; ++radius) {
            const u64 fast = solver.count_triangles_for_radius(radius);
            const u64 brute = brute_force_triangle_count_for_radius(radius);
            if (fast != brute) {
                std::cerr << "Validation failed for radius " << radius << ": fast=" << fast
                          << ", brute=" << brute << '\n';
                return false;
            }
        }
    }

    struct Checkpoint {
        int radius_limit;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {100, 4'950},
        {1200, 1'653'605},
    };

    for (const Checkpoint cp : checkpoints) {
        const Euler373Solver solver(cp.radius_limit);
        const u64 got = static_cast<u64>(solver.solve(false, 1));
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed for S(" << cp.radius_limit << "): got " << got
                      << ", expected " << cp.expected << '\n';
            return false;
        }
    }

    return true;
}

} // namespace

int main() {
    if (!run_validation_checkpoints()) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    const Euler373Solver solver(kTargetRadius);
    const u128 answer = solver.solve(true, threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
