#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 MOD = 101'001'001;
constexpr i64 BASE = 10;
constexpr u64 TARGET = 1'000'000'000'000'000'000ULL;

i64 normalize(i64 value) {
    value %= MOD;
    if (value < 0) {
        value += MOD;
    }
    return value;
}

i64 mod_pow(i64 base, u64 exponent) {
    i64 result = 1;
    base = normalize(base);
    while (exponent > 0) {
        if ((exponent & 1ULL) != 0) {
            result = result * base % MOD;
        }
        base = base * base % MOD;
        exponent >>= 1ULL;
    }
    return result;
}

i64 extended_gcd(const i64 a, const i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 next_x = 0;
    i64 next_y = 0;
    const i64 gcd = extended_gcd(b, a % b, next_x, next_y);
    x = next_y;
    y = next_x - (a / b) * next_y;
    return gcd;
}

i64 mod_inverse(const i64 value) {
    i64 x = 0;
    i64 y = 0;
    const i64 gcd = extended_gcd(normalize(value), MOD, x, y);
    if (gcd != 1) {
        std::cerr << "A required modular inverse does not exist.\n";
        std::exit(EXIT_FAILURE);
    }
    return normalize(x);
}

struct Node {
    u64 length = 0;
    int left = -1;
    int right = -1;
    int bit = -1;
};

struct Summary {
    i64 forward = 0;
    i64 reverse = 0;
    i64 pairs = 0;
    i64 ones = 0;
    bool ready = false;
};

struct QueryKey {
    int x = 0;
    int y = 0;
    i64 low = 0;
    i64 high = 0;
    int base_id = 0;

    bool operator==(const QueryKey& other) const {
        return x == other.x && y == other.y && low == other.low && high == other.high &&
               base_id == other.base_id;
    }
};

struct QueryKeyHash {
    std::size_t operator()(const QueryKey& key) const {
        std::size_t hash = static_cast<std::size_t>(key.x) * 1'000'003ULL +
                           static_cast<std::size_t>(key.y);
        hash ^= static_cast<std::size_t>(key.low) + 0x9e3779b97f4a7c15ULL + (hash << 6U) +
                (hash >> 2U);
        hash ^= static_cast<std::size_t>(key.high) + 0x9e3779b97f4a7c15ULL + (hash << 6U) +
                (hash >> 2U);
        hash ^= static_cast<std::size_t>(key.base_id) + (hash << 6U) + (hash >> 2U);
        return hash;
    }
};

class FibonacciSubwords {
  public:
    FibonacciSubwords() {
        bases_[0] = BASE;
        bases_[1] = mod_inverse(BASE);

        zero_ = add_bit(0);
        one_ = add_bit(1);
        fibonacci_.push_back(1);
        fibonacci_.push_back(2);
        standard_.push_back(zero_);
        standard_.push_back(concatenate(zero_, one_));

        while (fibonacci_.back() <= TARGET) {
            fibonacci_.push_back(fibonacci_[fibonacci_.size() - 1] +
                                 fibonacci_[fibonacci_.size() - 2]);
            standard_.push_back(concatenate(standard_[standard_.size() - 1],
                                            standard_[standard_.size() - 2]));
        }
    }

    i64 solve(const u64 k) {
        const auto iterator = std::upper_bound(fibonacci_.begin(), fibonacci_.end(), k);
        assert(iterator != fibonacci_.end());
        const int index = static_cast<int>(iterator - fibonacci_.begin());
        const u64 cycle_length = fibonacci_[static_cast<std::size_t>(index)];
        const u64 duplicate_count = cycle_length - 1 - k;
        const u64 gap_length = duplicate_count + 1;
        const int word = standard_[static_cast<std::size_t>(index)];

        const Summary word_base = summarize(word, 0);
        const Summary word_inverse = summarize(word, 1);
        const i64 cycle_power = power(0, cycle_length);
        const i64 inverse_cycle_power = power(1, cycle_length);

        const i64 all_correlations_base =
            normalize(word_base.pairs + cycle_power * word_inverse.pairs);
        const i64 all_correlations_inverse =
            normalize(word_inverse.pairs + inverse_cycle_power * word_base.pairs);

        const i64 prefix_correlations_base = normalize(
            difference_query(word, word, 1, static_cast<i64>(gap_length), 0) +
            cycle_power * difference_query(word,
                                           word,
                                           static_cast<i64>(cycle_length - gap_length),
                                           static_cast<i64>(cycle_length - 1),
                                           1));
        const i64 prefix_correlations_inverse = normalize(
            difference_query(word, word, 1, static_cast<i64>(gap_length), 1) +
            inverse_cycle_power * difference_query(word,
                                                   word,
                                                   static_cast<i64>(cycle_length - gap_length),
                                                   static_cast<i64>(cycle_length - 1),
                                                   0));

        const i64 used_correlations_base = normalize(
            all_correlations_base - cycle_power * prefix_correlations_inverse);
        const i64 used_correlations_inverse = normalize(
            all_correlations_inverse - inverse_cycle_power * prefix_correlations_base);

        const i64 inverse_geometric_denominator = mod_inverse(BASE * BASE - 1);
        const i64 geometric_squares =
            normalize(power(0, 2 * k) - 1) * inverse_geometric_denominator % MOD;
        const i64 full_cycle = normalize(
            word_base.ones * geometric_squares +
            2 * inverse_geometric_denominator % MOD *
                normalize(power(0, 2 * k) * used_correlations_inverse -
                          used_correlations_base));

        return normalize(full_cycle -
                         duplicated_window_sum(index, k, duplicate_count));
    }

  private:
    std::vector<Node> nodes_;
    std::vector<std::array<Summary, 2>> summary_cache_;
    std::map<std::pair<int, int>, int> concatenation_cache_;
    std::map<std::pair<int, u64>, int> prefix_cache_;
    std::vector<u64> fibonacci_;
    std::vector<int> standard_;
    std::array<i64, 2> bases_{};
    std::array<std::unordered_map<u64, i64>, 2> power_cache_;
    std::unordered_map<QueryKey, i64, QueryKeyHash> sum_query_cache_;
    std::unordered_map<QueryKey, i64, QueryKeyHash> difference_query_cache_;
    int zero_ = -1;
    int one_ = -1;

    int add_bit(const int bit) {
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{1, -1, -1, bit});
        summary_cache_.push_back({});
        return id;
    }

    int concatenate(const int left, const int right) {
        if (left < 0) {
            return right;
        }
        if (right < 0) {
            return left;
        }
        const std::pair<int, int> key{left, right};
        const auto found = concatenation_cache_.find(key);
        if (found != concatenation_cache_.end()) {
            return found->second;
        }
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{nodes_[static_cast<std::size_t>(left)].length +
                                  nodes_[static_cast<std::size_t>(right)].length,
                              left,
                              right,
                              -1});
        summary_cache_.push_back({});
        concatenation_cache_[key] = id;
        return id;
    }

    int prefix(const int standard_index, const u64 length) {
        assert(length > 0 && length <= fibonacci_[static_cast<std::size_t>(standard_index)]);
        if (length == fibonacci_[static_cast<std::size_t>(standard_index)]) {
            return standard_[static_cast<std::size_t>(standard_index)];
        }
        const std::pair<int, u64> key{standard_index, length};
        const auto found = prefix_cache_.find(key);
        if (found != prefix_cache_.end()) {
            return found->second;
        }
        assert(standard_index > 0);
        int result = -1;
        if (length <= fibonacci_[static_cast<std::size_t>(standard_index - 1)]) {
            result = prefix(standard_index - 1, length);
        } else {
            result = concatenate(
                standard_[static_cast<std::size_t>(standard_index - 1)],
                prefix(standard_index - 2,
                       length - fibonacci_[static_cast<std::size_t>(standard_index - 1)]));
        }
        prefix_cache_[key] = result;
        return result;
    }

    i64 power(const int base_id, const u64 exponent) {
        auto& cache = power_cache_[static_cast<std::size_t>(base_id)];
        const auto found = cache.find(exponent);
        if (found != cache.end()) {
            return found->second;
        }
        const i64 result = mod_pow(bases_[static_cast<std::size_t>(base_id)], exponent);
        cache[exponent] = result;
        return result;
    }

    i64 signed_power(const int base_id, const i64 exponent) {
        if (exponent >= 0) {
            return power(base_id, static_cast<u64>(exponent));
        }
        return power(1 - base_id, static_cast<u64>(-exponent));
    }

    Summary summarize(const int node_id, const int base_id) {
        Summary& cached =
            summary_cache_[static_cast<std::size_t>(node_id)][static_cast<std::size_t>(base_id)];
        if (cached.ready) {
            return cached;
        }

        const Node node = nodes_[static_cast<std::size_t>(node_id)];
        Summary result;
        if (node.bit >= 0) {
            result.forward = node.bit;
            result.reverse = node.bit;
            result.ones = node.bit;
            result.ready = true;
            cached = result;
            return result;
        }

        const Summary left = summarize(node.left, base_id);
        const Summary right = summarize(node.right, base_id);
        const u64 left_length = nodes_[static_cast<std::size_t>(node.left)].length;
        const u64 right_length = nodes_[static_cast<std::size_t>(node.right)].length;
        const i64 base = bases_[static_cast<std::size_t>(base_id)];

        result.forward = normalize(left.forward + power(base_id, left_length) * right.forward);
        result.reverse = normalize(power(base_id, right_length) * left.reverse + right.reverse);
        result.pairs =
            normalize(left.pairs + right.pairs + base * left.reverse % MOD * right.forward);
        result.ones = normalize(left.ones + right.ones);
        result.ready = true;
        cached = result;
        return result;
    }

    i64 sum_query(int x, int y, const i64 low, const i64 high, const int base_id) {
        if (x > y) {
            std::swap(x, y);
        }
        const QueryKey key{x, y, low, high, base_id};
        const auto found = sum_query_cache_.find(key);
        if (found != sum_query_cache_.end()) {
            return found->second;
        }

        const Node node_x = nodes_[static_cast<std::size_t>(x)];
        const Node node_y = nodes_[static_cast<std::size_t>(y)];
        const i64 maximum = static_cast<i64>(node_x.length + node_y.length - 2);
        i64 result = 0;

        if (high < 0 || low > maximum) {
            result = 0;
        } else if (low <= 0 && maximum <= high) {
            result = summarize(x, base_id).forward * summarize(y, base_id).forward % MOD;
        } else if (node_x.bit >= 0 && node_y.bit >= 0) {
            result = low <= 0 && 0 <= high ? node_x.bit * node_y.bit : 0;
        } else if (node_x.length >= node_y.length && node_x.bit < 0) {
            const u64 offset = nodes_[static_cast<std::size_t>(node_x.left)].length;
            result = normalize(
                sum_query(node_x.left, y, low, high, base_id) +
                power(base_id, offset) *
                    sum_query(node_x.right,
                              y,
                              low - static_cast<i64>(offset),
                              high - static_cast<i64>(offset),
                              base_id));
        } else {
            const u64 offset = nodes_[static_cast<std::size_t>(node_y.left)].length;
            result = normalize(
                sum_query(x, node_y.left, low, high, base_id) +
                power(base_id, offset) *
                    sum_query(x,
                              node_y.right,
                              low - static_cast<i64>(offset),
                              high - static_cast<i64>(offset),
                              base_id));
        }

        sum_query_cache_[key] = result;
        return result;
    }

    i64 difference_query(const int x,
                         const int y,
                         const i64 low,
                         const i64 high,
                         const int base_id) {
        const QueryKey key{x, y, low, high, base_id};
        const auto found = difference_query_cache_.find(key);
        if (found != difference_query_cache_.end()) {
            return found->second;
        }

        const Node node_x = nodes_[static_cast<std::size_t>(x)];
        const Node node_y = nodes_[static_cast<std::size_t>(y)];
        const i64 minimum = -static_cast<i64>(node_x.length - 1);
        const i64 maximum = static_cast<i64>(node_y.length - 1);
        i64 result = 0;

        if (high < minimum || low > maximum) {
            result = 0;
        } else if (low <= minimum && maximum <= high) {
            result = summarize(x, 1 - base_id).forward * summarize(y, base_id).forward % MOD;
        } else if (node_x.bit >= 0 && node_y.bit >= 0) {
            result = low <= 0 && 0 <= high ? node_x.bit * node_y.bit : 0;
        } else if (node_x.length >= node_y.length && node_x.bit < 0) {
            const u64 offset = nodes_[static_cast<std::size_t>(node_x.left)].length;
            result = normalize(
                difference_query(node_x.left, y, low, high, base_id) +
                signed_power(base_id, -static_cast<i64>(offset)) *
                    difference_query(node_x.right,
                                     y,
                                     low + static_cast<i64>(offset),
                                     high + static_cast<i64>(offset),
                                     base_id));
        } else {
            const u64 offset = nodes_[static_cast<std::size_t>(node_y.left)].length;
            result = normalize(
                difference_query(x, node_y.left, low, high, base_id) +
                power(base_id, offset) *
                    difference_query(x,
                                     node_y.right,
                                     low - static_cast<i64>(offset),
                                     high - static_cast<i64>(offset),
                                     base_id));
        }

        difference_query_cache_[key] = result;
        return result;
    }

    i64 duplicated_window_sum(const int standard_index,
                              const u64 k,
                              const u64 count) {
        if (count == 0) {
            return 0;
        }

        const int initial_word = prefix(standard_index, k);
        const i64 initial_value = summarize(initial_word, 0).reverse;
        if (count == 1) {
            return initial_value * initial_value % MOD;
        }

        const u64 transition_count = count - 1;
        const int outgoing = prefix(standard_index, transition_count);
        const Summary outgoing_summary = summarize(outgoing, 0);
        const i64 window_power = power(0, k);
        const i64 forward_delta = normalize(outgoing_summary.reverse -
                                            window_power * outgoing_summary.forward);
        const i64 reverse_delta = normalize(outgoing_summary.forward -
                                            window_power * outgoing_summary.reverse);

        const i64 low_inverse = sum_query(outgoing,
                                          outgoing,
                                          0,
                                          static_cast<i64>(transition_count) - 2,
                                          1);
        const i64 high_base = sum_query(outgoing,
                                        outgoing,
                                        static_cast<i64>(transition_count),
                                        static_cast<i64>(2 * transition_count - 2),
                                        0);
        const i64 diagonal_base = sum_query(outgoing,
                                            outgoing,
                                            static_cast<i64>(transition_count - 1),
                                            static_cast<i64>(transition_count - 1),
                                            0);

        const i64 diagonal =
            diagonal_base * power(1, transition_count - 1) % MOD;
        const i64 lower_cross = power(0, transition_count - 1) * low_inverse % MOD;
        const i64 upper_cross = power(1, transition_count - 1) * high_base % MOD;
        const i64 delta_pairs = normalize(
            normalize(1 + window_power * window_power) * outgoing_summary.pairs -
            window_power * normalize(lower_cross + upper_cross));
        const i64 delta_squares = normalize(
            normalize(1 + window_power * window_power) * outgoing_summary.ones -
            2 * window_power % MOD * diagonal);

        const i64 last_value = normalize(
            power(0, transition_count) * initial_value + reverse_delta);
        const i64 value_delta_sum = normalize(
            initial_value * forward_delta + bases_[1] * delta_pairs);
        const i64 numerator = normalize(
            initial_value * initial_value -
            BASE * BASE % MOD * last_value % MOD * last_value +
            2 * BASE % MOD * value_delta_sum + delta_squares);

        return numerator * mod_inverse(1 - BASE * BASE) % MOD;
    }
};

i64 brute_force(const int k) {
    std::string older = "0";
    std::string newer = "01";
    while (static_cast<int>(newer.size()) < 10 * k + 20) {
        const std::string next = newer + older;
        older = newer;
        newer = next;
    }

    std::set<std::string> factors;
    for (int start = 0; start + k <= static_cast<int>(newer.size()); ++start) {
        factors.insert(newer.substr(static_cast<std::size_t>(start), static_cast<std::size_t>(k)));
    }
    assert(factors.size() == static_cast<std::size_t>(k + 1));

    i64 result = 0;
    for (const std::string& factor : factors) {
        i64 value = 0;
        for (const char digit : factor) {
            value = (BASE * value + digit - '0') % MOD;
        }
        result = (result + value * value) % MOD;
    }
    return result;
}

void run_checkpoints() {
    FibonacciSubwords solver;
    if (brute_force(3) != 20'302) {
        std::cerr << "Checkpoint failed: Psi(3).\n";
        std::exit(EXIT_FAILURE);
    }
    for (int k = 1; k <= 50; ++k) {
        if (solver.solve(static_cast<u64>(k)) != brute_force(k)) {
            std::cerr << "Checkpoint failed: brute force comparison for k=" << k << ".\n";
            std::exit(EXIT_FAILURE);
        }
    }
    if (solver.solve(10) != 10'699'667) {
        std::cerr << "Checkpoint failed: supplied Psi(10).\n";
        std::exit(EXIT_FAILURE);
    }
    std::cerr << "Validation checkpoints passed.\n";
}

}

int main() {
    run_checkpoints();
    FibonacciSubwords solver;
    std::cout << solver.solve(TARGET) << '\n';
    return 0;
}
