#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kMod = 1'000'000'993U;
constexpr u32 kDefaultN = 11'111'111U;

struct Options {
    u32 n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

inline u32 mul_mod(u32 a, u32 b) {
    return static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(b)) % kMod);
}

inline u32 add_mod(u32 a, u32 b) {
    const u32 sum = a + b;
    return (sum >= kMod) ? (sum - kMod) : sum;
}

u32 mod_pow(u32 base, u32 exponent) {
    u32 result = 1U;
    u32 b = base;
    u32 e = exponent;
    while (e > 0U) {
        if ((e & 1U) != 0U) {
            result = mul_mod(result, b);
        }
        b = mul_mod(b, b);
        e >>= 1U;
    }
    return result;
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
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
        const u32 digit = static_cast<u32>(c - '0');
        parsed = parsed * 10ULL + static_cast<u64>(digit);
        if (parsed > static_cast<u64>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
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

        u32 parsed_u32 = 0;
        if (parse_u32_after_prefix(arg, "--n=", parsed_u32)) {
            options.n = parsed_u32;
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n == 0U) {
        std::cerr << "--n must be >= 1.\n";
        return false;
    }

    return true;
}

u64 largest_b_smooth_divisor(u32 b, u64 x) {
    u64 result = 1ULL;
    u64 value = x;

    for (u64 p = 2ULL; p * p <= value; ++p) {
        if (value % p != 0ULL) {
            continue;
        }

        u64 prime_power = 1ULL;
        while (value % p == 0ULL) {
            value /= p;
            prime_power *= p;
        }

        if (p <= static_cast<u64>(b)) {
            result *= prime_power;
        }
    }

    if (value > 1ULL && value <= static_cast<u64>(b)) {
        result *= value;
    }

    return result;
}

struct SieveData {
    u32 n = 0U;
    std::vector<u32> weights;
    std::vector<u32> strip_next;
    std::vector<u32> strip_pos;
    std::vector<u32> strip_mul;
    std::vector<u32> strip_inv;
};

SieveData build_sieve_data(u32 n) {
    SieveData data;
    data.n = n;
    std::vector<u32> spf(static_cast<std::size_t>(n) + 1ULL, 0U);
    std::vector<int> prime_pos(static_cast<std::size_t>(n) + 1ULL, -1);
    std::vector<u32> primes;

    primes.reserve(static_cast<std::size_t>(n / 10U));

    for (u32 i = 2U; i <= n; ++i) {
        if (spf[i] == 0U) {
            spf[i] = i;
            primes.push_back(i);
        }

        for (const u32 p : primes) {
            const u64 composite = static_cast<u64>(i) * static_cast<u64>(p);
            if (composite > static_cast<u64>(n)) {
                break;
            }
            spf[static_cast<std::size_t>(composite)] = p;
            if (p == spf[i]) {
                break;
            }
        }
    }

    std::vector<u32> inv_prime_by_pos(primes.size());
    data.weights.resize(primes.size());
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const u32 p = primes[i];
        prime_pos[p] = static_cast<int>(i);
        inv_prime_by_pos[i] = mod_pow(p, kMod - 2U);
        const u32 next_prime = (i + 1U < primes.size()) ? primes[i + 1U] : (n + 1U);
        data.weights[i] = next_prime - p;
    }

    data.strip_next.assign(static_cast<std::size_t>(n) + 1ULL, 1U);
    data.strip_pos.assign(static_cast<std::size_t>(n) + 1ULL, 0U);
    data.strip_mul.assign(static_cast<std::size_t>(n) + 1ULL, 1U);
    data.strip_inv.assign(static_cast<std::size_t>(n) + 1ULL, 1U);

    for (u32 x = 2U; x <= n; ++x) {
        const u32 p = spf[x];
        const u32 pos = static_cast<u32>(prime_pos[p]);
        const u32 inv_p = inv_prime_by_pos[static_cast<std::size_t>(pos)];

        u32 value = x;
        u32 mul = 1U;
        u32 inv_mul = 1U;
        do {
            value /= p;
            mul = mul_mod(mul, p);
            inv_mul = mul_mod(inv_mul, inv_p);
        } while (value % p == 0U);

        data.strip_next[x] = value;
        data.strip_pos[x] = pos;
        data.strip_mul[x] = mul;
        data.strip_inv[x] = inv_mul;
    }

    return data;
}

class SegmentTree {
  public:
    explicit SegmentTree(const std::vector<u32>& weights) {
        leaf_count_ = static_cast<int>(weights.size());
        size_ = 1;
        while (size_ < leaf_count_) {
            size_ <<= 1;
        }

        nodes_.assign(static_cast<std::size_t>(2 * size_), Node{1U, 0U});

        for (int i = 0; i < leaf_count_; ++i) {
            nodes_[static_cast<std::size_t>(size_ + i)].sum =
                weights[static_cast<std::size_t>(i)] % kMod;
        }

        for (int node = size_ - 1; node >= 1; --node) {
            pull(node);
        }
    }

    void multiply_leaf(int leaf_index, u32 factor) {
        int node = size_ + leaf_index;
        Node& leaf = nodes_[static_cast<std::size_t>(node)];
        leaf.prod = mul_mod(leaf.prod, factor);
        leaf.sum = mul_mod(leaf.sum, factor);

        for (node >>= 1; node > 0; node >>= 1) {
            pull(node);
        }
    }

    u32 root_sum() const {
        return nodes_[1].sum;
    }

  private:
    struct Node {
        u32 prod;
        u32 sum;
    };

    void pull(int node) {
        const int left = node << 1;
        const int right = left | 1;
        const Node& lhs = nodes_[static_cast<std::size_t>(left)];
        const Node& rhs = nodes_[static_cast<std::size_t>(right)];
        Node& out = nodes_[static_cast<std::size_t>(node)];
        out.prod = mul_mod(lhs.prod, rhs.prod);
        out.sum = add_mod(lhs.sum, mul_mod(lhs.prod, rhs.sum));
    }

    int leaf_count_ = 0;
    int size_ = 0;
    std::vector<Node> nodes_;
};

class SolverState {
  public:
    explicit SolverState(const SieveData& data) : data_(data), tree_(data.weights) {}

    u32 current_row_contribution() const {
        const u32 v = tree_.root_sum();
        return (v + 1U == kMod) ? 0U : (v + 1U);
    }

    void apply_ratio(u32 numerator, u32 denominator) {
        u32 positions[16];
        u32 factors[16];
        int used = 0;

        auto append = [&](u32 value, bool divide) {
            while (value > 1U) {
                const u32 pos = data_.strip_pos[value];
                const u32 factor = divide ? data_.strip_inv[value] : data_.strip_mul[value];

                int slot = -1;
                for (int i = 0; i < used; ++i) {
                    if (positions[i] == pos) {
                        slot = i;
                        break;
                    }
                }

                if (slot < 0) {
                    positions[used] = pos;
                    factors[used] = factor;
                    ++used;
                } else {
                    factors[slot] = mul_mod(factors[slot], factor);
                }
                value = data_.strip_next[value];
            }
        };

        append(numerator, false);
        append(denominator, true);

        for (int i = 0; i < used; ++i) {
            if (factors[i] != 1U) {
                tree_.multiply_leaf(static_cast<int>(positions[i]), factors[i]);
            }
        }
    }

  private:
    const SieveData& data_;
    SegmentTree tree_;
};

u32 sum_symmetric(const SieveData& data, u32 n) {
    SolverState state(data);
    u32 total = 0U;
    const u32 half = n / 2U;

    for (u32 r = 0U;; ++r) {
        const u32 contribution = state.current_row_contribution();
        if (r == n - r) {
            total = add_mod(total, contribution);
        } else {
            total = add_mod(total, add_mod(contribution, contribution));
        }

        if (r == half) {
            break;
        }
        state.apply_ratio(n - r, r + 1U);
    }

    return total;
}

u32 solve_mod(const SieveData& data,
              u32 n,
              bool allow_multithreading,
              unsigned requested_threads) {
    (void)allow_multithreading;
    (void)requested_threads;
    return sum_symmetric(data, n);
}

u32 solve_for_n(u32 n, bool allow_multithreading, unsigned requested_threads) {
    SieveData data = build_sieve_data(n);
    return solve_mod(data, n, allow_multithreading, requested_threads);
}

bool expect_equal_u64(const char* label, u64 actual, u64 expected) {
    if (actual == expected) {
        return true;
    }

    std::cerr << "Checkpoint failed for " << label << ": got " << actual << ", expected "
              << expected << '\n';
    return false;
}

bool expect_equal_u32(const char* label, u32 actual, u32 expected) {
    return expect_equal_u64(label, static_cast<u64>(actual), static_cast<u64>(expected));
}

bool run_checkpoints() {
    bool ok = true;

    ok = ok && expect_equal_u64("S_1(10)", largest_b_smooth_divisor(1U, 10ULL), 1ULL);
    ok = ok && expect_equal_u64("S_4(2100)", largest_b_smooth_divisor(4U, 2100ULL), 12ULL);
    ok = ok && expect_equal_u64("S_17(2496144)", largest_b_smooth_divisor(17U, 2'496'144ULL),
                                5'712ULL);

    ok = ok && expect_equal_u32("F(11)", solve_for_n(11U, false, 1U), 3'132U);
    ok = ok &&
         expect_equal_u32("F(1111) mod M", solve_for_n(1'111U, false, 1U), 706'036'312U);
    ok = ok &&
         expect_equal_u32("F(111111) mod M", solve_for_n(111'111U, false, 1U), 22'156'169U);

    if (!ok) {
        std::cerr << "At least one checkpoint failed.\n";
    }
    return ok;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }

    const auto start = std::chrono::steady_clock::now();
    const SieveData data = build_sieve_data(options.n);
    const u32 answer =
        solve_mod(data, options.n, options.allow_multithreading, options.requested_threads);
    const auto end = std::chrono::steady_clock::now();

    const std::chrono::duration<long double> elapsed = end - start;
    std::cout << "F(" << options.n << ") mod " << kMod << " = " << answer << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " seconds\n";

    return 0;
}
