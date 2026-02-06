#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultOrder = 35;

struct Options {
    u64 order = kDefaultOrder;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct Fraction {
    u64 num = 0;
    u64 den = 1;

    bool operator==(const Fraction& other) const noexcept {
        return num == other.num && den == other.den;
    }
};

struct FractionHash {
    std::size_t operator()(const Fraction& f) const noexcept {
        const u64 x = f.num * 0x9E3779B185EBCA87ULL;
        const u64 y = f.den + 0xC2B2AE3D27D4EB4FULL;
        return static_cast<std::size_t>(x ^ (y + (x << 6) + (x >> 2)));
    }
};

u128 gcd_u128(u128 a, u128 b) {
    while (b != 0) {
        const u128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

u64 cast_u64_checked(u128 value) {
    if (value > static_cast<u128>(std::numeric_limits<u64>::max())) {
        std::cerr << "Overflow while converting to uint64_t.\n";
        std::exit(1);
    }
    return static_cast<u64>(value);
}

Fraction make_fraction(u64 num, u64 den) {
    if (den == 0) {
        std::cerr << "Attempted to build fraction with zero denominator.\n";
        std::exit(1);
    }
    const u64 g = std::gcd(num, den);
    return Fraction{num / g, den / g};
}

Fraction add_fraction(const Fraction& a, const Fraction& b) {
    const u64 g = std::gcd(a.den, b.den);
    u128 num = static_cast<u128>(a.num) * static_cast<u128>(b.den / g) +
               static_cast<u128>(b.num) * static_cast<u128>(a.den / g);
    u128 den = static_cast<u128>(a.den / g) * static_cast<u128>(b.den);
    const u128 h = gcd_u128(num, den);
    num /= h;
    den /= h;
    return Fraction{cast_u64_checked(num), cast_u64_checked(den)};
}

Fraction mul_fraction(const Fraction& a, const Fraction& b) {
    u128 num = static_cast<u128>(a.num) * static_cast<u128>(b.num);
    u128 den = static_cast<u128>(a.den) * static_cast<u128>(b.den);
    const u128 g = gcd_u128(num, den);
    num /= g;
    den /= g;
    return Fraction{cast_u64_checked(num), cast_u64_checked(den)};
}

Fraction reciprocal(const Fraction& f) {
    if (f.num == 0) {
        std::cerr << "Attempted reciprocal of zero.\n";
        std::exit(1);
    }
    return Fraction{f.den, f.num};
}

Fraction div_fraction(const Fraction& a, const Fraction& b) {
    return mul_fraction(a, reciprocal(b));
}

Fraction square_fraction(const Fraction& f) {
    return mul_fraction(f, f);
}

bool is_square_u64(u64 value, u64& root) {
    const long double approx = std::sqrt(static_cast<long double>(value));
    u64 r = static_cast<u64>(approx);
    while (static_cast<u128>(r) * static_cast<u128>(r) < value) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > value) {
        --r;
    }
    if (static_cast<u128>(r) * static_cast<u128>(r) == value) {
        root = r;
        return true;
    }
    return false;
}

bool sqrt_fraction(const Fraction& f, Fraction& root) {
    u64 nroot = 0;
    u64 droot = 0;
    if (!is_square_u64(f.num, nroot) || !is_square_u64(f.den, droot)) {
        return false;
    }
    root = make_fraction(nroot, droot);
    return true;
}

std::vector<Fraction> build_fraction_list(u64 order) {
    std::vector<Fraction> fractions;
    for (u64 den = 2; den <= order; ++den) {
        for (u64 num = 1; num < den; ++num) {
            if (std::gcd(num, den) == 1) {
                fractions.push_back(Fraction{num, den});
            }
        }
    }
    return fractions;
}

unsigned decide_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t x_count) {
    if (!allow_multithreading) {
        return 1U;
    }
    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }
    if (x_count < 128U) {
        threads = 1U;
    }
    if (threads > x_count) {
        threads = static_cast<unsigned>(x_count);
    }
    return std::max(1U, threads);
}

void add_sum_if_member(const Fraction& x,
                       const Fraction& y,
                       const Fraction& z,
                       const std::unordered_set<Fraction, FractionHash>& allowed,
                       std::unordered_set<Fraction, FractionHash>& sums) {
    if (allowed.find(z) == allowed.end()) {
        return;
    }
    const Fraction s = add_fraction(add_fraction(x, y), z);
    sums.insert(s);
}

void process_pair(const Fraction& x,
                  const Fraction& y,
                  const std::unordered_set<Fraction, FractionHash>& allowed,
                  std::unordered_set<Fraction, FractionHash>& sums) {
    // For positive rationals, f_n(x,y,z)=0 is equivalent to x^n + y^n = z^n.
    // Non-trivial rational solutions only occur for n in {-2,-1,1,2}.

    // n = 1: x + y = z.
    const Fraction z1 = add_fraction(x, y);
    add_sum_if_member(x, y, z1, allowed, sums);

    // n = -1: 1/x + 1/y = 1/z  =>  z = xy/(x+y).
    const Fraction zminus1 = div_fraction(mul_fraction(x, y), z1);
    add_sum_if_member(x, y, zminus1, allowed, sums);

    // n = 2: x^2 + y^2 = z^2.
    const Fraction z2_sq = add_fraction(square_fraction(x), square_fraction(y));
    Fraction z2{};
    if (sqrt_fraction(z2_sq, z2)) {
        add_sum_if_member(x, y, z2, allowed, sums);
    }

    // n = -2: 1/x^2 + 1/y^2 = 1/z^2  =>  z^2 = 1 / (1/x^2 + 1/y^2).
    const Fraction inv_x_sq = square_fraction(reciprocal(x));
    const Fraction inv_y_sq = square_fraction(reciprocal(y));
    const Fraction zminus2_sq = reciprocal(add_fraction(inv_x_sq, inv_y_sq));
    Fraction zminus2{};
    if (sqrt_fraction(zminus2_sq, zminus2)) {
        add_sum_if_member(x, y, zminus2, allowed, sums);
    }
}

std::unordered_set<Fraction, FractionHash> collect_sums_fast(
    u64 order,
    bool allow_multithreading,
    unsigned requested_threads) {
    const std::vector<Fraction> fractions = build_fraction_list(order);
    const std::unordered_set<Fraction, FractionHash> allowed(fractions.begin(), fractions.end());
    const std::size_t n = fractions.size();

    const unsigned threads = decide_thread_count(allow_multithreading, requested_threads, n);
    if (threads <= 1U) {
        std::unordered_set<Fraction, FractionHash> sums;
        sums.reserve(n * 4);
        for (const Fraction& x : fractions) {
            for (const Fraction& y : fractions) {
                process_pair(x, y, allowed, sums);
            }
        }
        return sums;
    }

    std::vector<std::unordered_set<Fraction, FractionHash>> locals(threads);
    for (auto& local : locals) {
        local.reserve((n * 4) / threads + 64U);
    }

    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        const std::size_t begin = (n * tid) / threads;
        const std::size_t end = (n * (tid + 1U)) / threads;

        workers.emplace_back([&, tid, begin, end]() {
            auto& local = locals[tid];
            for (std::size_t i = begin; i < end; ++i) {
                const Fraction& x = fractions[i];
                for (const Fraction& y : fractions) {
                    process_pair(x, y, allowed, local);
                }
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    std::unordered_set<Fraction, FractionHash> sums;
    for (const auto& local : locals) {
        sums.insert(local.begin(), local.end());
    }
    return sums;
}

bool satisfies_order_from_derived_equations(const Fraction& x,
                                            const Fraction& y,
                                            const Fraction& z) {
    if (add_fraction(x, y) == z) {
        return true;
    }

    if (div_fraction(mul_fraction(x, y), add_fraction(x, y)) == z) {
        return true;
    }

    if (add_fraction(square_fraction(x), square_fraction(y)) == square_fraction(z)) {
        return true;
    }

    const Fraction inv_x_sq = square_fraction(reciprocal(x));
    const Fraction inv_y_sq = square_fraction(reciprocal(y));
    const Fraction inv_z_sq = square_fraction(reciprocal(z));
    if (add_fraction(inv_x_sq, inv_y_sq) == inv_z_sq) {
        return true;
    }

    return false;
}

struct RationalBig {
    cpp_int num = 0;
    cpp_int den = 1;

    RationalBig() = default;
    RationalBig(cpp_int n, cpp_int d) : num(std::move(n)), den(std::move(d)) {
        normalize();
    }

    static cpp_int abs_value(cpp_int v) {
        return (v < 0) ? -v : v;
    }

    static cpp_int gcd(cpp_int a, cpp_int b) {
        a = abs_value(a);
        b = abs_value(b);
        while (b != 0) {
            cpp_int r = a % b;
            a = b;
            b = r;
        }
        return a;
    }

    void normalize() {
        if (den == 0) {
            std::cerr << "Big rational with zero denominator.\n";
            std::exit(1);
        }
        if (den < 0) {
            den = -den;
            num = -num;
        }
        const cpp_int g = gcd(num, den);
        num /= g;
        den /= g;
    }
};

RationalBig operator+(const RationalBig& a, const RationalBig& b) {
    return RationalBig(a.num * b.den + b.num * a.den, a.den * b.den);
}

RationalBig operator-(const RationalBig& a, const RationalBig& b) {
    return RationalBig(a.num * b.den - b.num * a.den, a.den * b.den);
}

RationalBig operator*(const RationalBig& a, const RationalBig& b) {
    return RationalBig(a.num * b.num, a.den * b.den);
}

RationalBig operator/(const RationalBig& a, const RationalBig& b) {
    return RationalBig(a.num * b.den, a.den * b.num);
}

RationalBig pow_rational(RationalBig base, int exponent) {
    if (exponent == 0) {
        return RationalBig(1, 1);
    }
    if (exponent < 0) {
        base = RationalBig(base.den, base.num);
        exponent = -exponent;
    }

    RationalBig result(1, 1);
    int e = exponent;
    while (e > 0) {
        if ((e & 1) != 0) {
            result = result * base;
        }
        base = base * base;
        e >>= 1;
    }
    return result;
}

RationalBig to_big(const Fraction& f) {
    return RationalBig(cpp_int(f.num), cpp_int(f.den));
}

bool satisfies_order_from_original_formula(const Fraction& x,
                                           const Fraction& y,
                                           const Fraction& z) {
    static const int exponents[] = {-2, -1, 1, 2};

    const RationalBig bx = to_big(x);
    const RationalBig by = to_big(y);
    const RationalBig bz = to_big(z);

    for (const int n : exponents) {
        const RationalBig f1 = pow_rational(bx, n + 1) + pow_rational(by, n + 1) -
                               pow_rational(bz, n + 1);

        const RationalBig symmetric = bx * by + by * bz + bz * bx;
        const RationalBig f2 = symmetric * (pow_rational(bx, n - 1) +
                                            pow_rational(by, n - 1) -
                                            pow_rational(bz, n - 1));

        const RationalBig f3 = bx * by * bz * (pow_rational(bx, n - 2) +
                                               pow_rational(by, n - 2) -
                                               pow_rational(bz, n - 2));

        const RationalBig fn = f1 + f2 - f3;
        if (fn.num == 0) {
            return true;
        }
    }
    return false;
}

std::unordered_set<Fraction, FractionHash> collect_sums_bruteforce_derived(u64 order) {
    const std::vector<Fraction> fractions = build_fraction_list(order);
    std::unordered_set<Fraction, FractionHash> sums;
    sums.reserve(fractions.size() * 8);

    for (const Fraction& x : fractions) {
        for (const Fraction& y : fractions) {
            for (const Fraction& z : fractions) {
                if (satisfies_order_from_derived_equations(x, y, z)) {
                    sums.insert(add_fraction(add_fraction(x, y), z));
                }
            }
        }
    }
    return sums;
}

std::unordered_set<Fraction, FractionHash> collect_sums_bruteforce_original(u64 order) {
    const std::vector<Fraction> fractions = build_fraction_list(order);
    std::unordered_set<Fraction, FractionHash> sums;
    sums.reserve(fractions.size() * 8);

    for (const Fraction& x : fractions) {
        for (const Fraction& y : fractions) {
            for (const Fraction& z : fractions) {
                if (satisfies_order_from_original_formula(x, y, z)) {
                    sums.insert(add_fraction(add_fraction(x, y), z));
                }
            }
        }
    }
    return sums;
}

bool set_equal(const std::unordered_set<Fraction, FractionHash>& a,
               const std::unordered_set<Fraction, FractionHash>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (const Fraction& value : a) {
        if (b.find(value) == b.end()) {
            return false;
        }
    }
    return true;
}

std::pair<cpp_int, cpp_int> sum_distinct_s_values(
    const std::unordered_set<Fraction, FractionHash>& sums) {
    cpp_int num = 0;
    cpp_int den = 1;

    for (const Fraction& f : sums) {
        num = num * f.den + den * f.num;
        den *= f.den;
        const cpp_int g = RationalBig::gcd(num, den);
        num /= g;
        den /= g;
    }
    return {num, den};
}

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }
    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        u64 parsed_order = 0;
        if (parse_u64_after_prefix(arg, "--order=", parsed_order)) {
            options.order = parsed_order;
            continue;
        }

        unsigned parsed_threads = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_threads)) {
            options.requested_threads = parsed_threads;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.order < 2) {
        std::cerr << "--order must be at least 2.\n";
        return false;
    }

    return true;
}

bool run_validation_checkpoints(const Options& options) {
    // Checkpoint 1: direct original-formula brute force vs derived equations.
    {
        const u64 k = 7;
        const auto original = collect_sums_bruteforce_original(k);
        const auto derived = collect_sums_bruteforce_derived(k);
        if (!set_equal(original, derived)) {
            std::cerr << "Checkpoint failed: original formula and derived equations differ at "
                      << "order " << k << ".\n";
            return false;
        }
    }

    // Checkpoint 2: fast pair-scan vs brute-force triples.
    {
        const u64 k = 12;
        const auto brute = collect_sums_bruteforce_derived(k);
        const auto fast = collect_sums_fast(k, false, 1);
        if (!set_equal(brute, fast)) {
            std::cerr << "Checkpoint failed: fast solver mismatch at order " << k << ".\n";
            return false;
        }
    }

    // Checkpoint 3: single-thread vs multi-thread consistency on target order.
    {
        const auto single = collect_sums_fast(options.order, false, 1);
        const auto multi = collect_sums_fast(options.order,
                                             options.allow_multithreading,
                                             options.requested_threads);
        if (!set_equal(single, multi)) {
            std::cerr << "Checkpoint failed: thread-consistency mismatch at order "
                      << options.order << ".\n";
            return false;
        }
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
        if (!run_validation_checkpoints(options)) {
            return 1;
        }
    }

    const auto sums = collect_sums_fast(options.order,
                                        options.allow_multithreading,
                                        options.requested_threads);
    const auto [u, v] = sum_distinct_s_values(sums);
    std::cout << (u + v) << '\n';
    return 0;
}
