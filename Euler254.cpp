#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int max_i = 150;
    bool run_checkpoints = true;
};

constexpr int kMaxLen = 18;
constexpr int kMaxSum = 150;
constexpr int kMod = 362880;            // 9!
constexpr int kWordBits = 64;
constexpr int kWords = kMod / kWordBits;  // exact: 5670

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--max-i=", options.max_i)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_i >= 1 && options.max_i <= kMaxSum;
}

struct Node {
    std::array<u64, 10> cnt{};
    bool valid = false;
};

u64 digit_sum_of_node(const Node& node) {
    u64 sum = 0;
    for (int d = 1; d <= 9; ++d) {
        sum += static_cast<u64>(d) * node.cnt[static_cast<std::size_t>(d)];
    }
    return sum;
}

u64 length_of_node(const Node& node) {
    return std::accumulate(node.cnt.begin(), node.cnt.end(), 0ULL);
}

bool node_less(const Node& lhs, const Node& rhs) {
    if (!rhs.valid) {
        return lhs.valid;
    }
    if (!lhs.valid) {
        return false;
    }

    const u64 len_l = length_of_node(lhs);
    const u64 len_r = length_of_node(rhs);
    if (len_l != len_r) {
        return len_l < len_r;
    }

    for (int d = 0; d <= 9; ++d) {
        const u64 cl = lhs.cnt[static_cast<std::size_t>(d)];
        const u64 cr = rhs.cnt[static_cast<std::size_t>(d)];
        if (cl != cr) {
            return cl > cr;
        }
    }

    return false;
}

class DpTable {
public:
    DpTable() {
        data_.assign(static_cast<std::size_t>(kMaxLen + 1) * static_cast<std::size_t>(kMaxSum + 1) *
                         static_cast<std::size_t>(kWords),
                     0ULL);
        nonempty_.assign(static_cast<std::size_t>(kMaxLen + 1) * static_cast<std::size_t>(kMaxSum + 1), false);
    }

    u64* ptr(int len, int sum) {
        return data_.data() + offset(len, sum);
    }

    const u64* ptr(int len, int sum) const {
        return data_.data() + offset(len, sum);
    }

    bool nonempty(int len, int sum) const {
        return nonempty_[static_cast<std::size_t>(len) * static_cast<std::size_t>(kMaxSum + 1) +
                         static_cast<std::size_t>(sum)];
    }

    void set_nonempty(int len, int sum) {
        nonempty_[static_cast<std::size_t>(len) * static_cast<std::size_t>(kMaxSum + 1) +
                  static_cast<std::size_t>(sum)] = true;
    }

    static void set_bit(u64* bits, int pos) {
        bits[pos / 64] |= 1ULL << (pos % 64);
    }

    static bool get_bit(const u64* bits, int pos) {
        return ((bits[pos / 64] >> (pos % 64)) & 1ULL) != 0ULL;
    }

private:
    std::size_t offset(int len, int sum) const {
        return (static_cast<std::size_t>(len) * static_cast<std::size_t>(kMaxSum + 1) +
                static_cast<std::size_t>(sum)) *
               static_cast<std::size_t>(kWords);
    }

    std::vector<u64> data_;
    std::vector<bool> nonempty_;
};

void rotate_or(const u64* src, u64* dst, int shift) {
    if (shift < 0) {
        shift %= kMod;
        if (shift < 0) {
            shift += kMod;
        }
    } else if (shift >= kMod) {
        shift %= kMod;
    }

    if (shift == 0) {
        for (int i = 0; i < kWords; ++i) {
            dst[i] |= src[i];
        }
        return;
    }

    const int word_shift = shift / 64;
    const int bit_shift = shift % 64;

    if (bit_shift == 0) {
        for (int i = 0; i < kWords; ++i) {
            dst[(i + word_shift) % kWords] |= src[i];
        }
        return;
    }

    for (int i = 0; i < kWords; ++i) {
        const u64 v = src[i];
        if (v == 0) {
            continue;
        }
        const int j = (i + word_shift) % kWords;
        dst[j] |= v << bit_shift;
        dst[(j + 1) % kWords] |= v >> (64 - bit_shift);
    }
}

class Solver {
public:
    Solver() {
        pow10mod_.fill(0);
        pow10mod_[0] = 1;
        for (int i = 1; i <= kMaxLen; ++i) {
            pow10mod_[static_cast<std::size_t>(i)] =
                (pow10mod_[static_cast<std::size_t>(i - 1)] * 10) % kMod;
        }
        build_dp();
    }

    u64 solve_sum_sg(const int max_i) {
        u64 total = 0;
        for (int i = 1; i <= max_i; ++i) {
            Node best = best_node_for_sf(i);
            total += digit_sum_of_node(best);
        }
        return total;
    }

    Node best_node_for_sf(const int sf_value) {
        Node best;

        for (int residue = 0; residue < kMod; ++residue) {
            int min_len = -1;
            for (int len = 1; len <= kMaxLen; ++len) {
                if (!dp_.nonempty(len, sf_value)) {
                    continue;
                }
                if (DpTable::get_bit(dp_.ptr(len, sf_value), residue)) {
                    min_len = len;
                    break;
                }
            }
            if (min_len < 0) {
                continue;
            }

            const u64 y = backtrace_min_number(min_len, sf_value, residue);
            const Node candidate = node_from_y(y);
            if (node_less(candidate, best)) {
                best = candidate;
            }
        }

        best.valid = true;
        return best;
    }

    u64 backtrace_min_number(int len, int sum, int residue) const {
        u64 value = 0;
        while (len > 0) {
            bool found = false;
            for (int d = 0; d <= 9 && d <= sum; ++d) {
                int prev_residue = residue - (d * pow10mod_[static_cast<std::size_t>(len - 1)]) % kMod;
                prev_residue %= kMod;
                if (prev_residue < 0) {
                    prev_residue += kMod;
                }

                if (!dp_.nonempty(len - 1, sum - d)) {
                    continue;
                }
                if (!DpTable::get_bit(dp_.ptr(len - 1, sum - d), prev_residue)) {
                    continue;
                }

                value = value * 10ULL + static_cast<u64>(d);
                sum -= d;
                residue = prev_residue;
                --len;
                found = true;
                break;
            }

            if (!found) {
                return std::numeric_limits<u64>::max();
            }
        }
        return value;
    }

private:
    void build_dp() {
        DpTable::set_bit(dp_.ptr(0, 0), 0);
        dp_.set_nonempty(0, 0);

        for (int len = 0; len < kMaxLen; ++len) {
            for (int sum = 0; sum <= kMaxSum; ++sum) {
                if (!dp_.nonempty(len, sum)) {
                    continue;
                }
                const u64* src = dp_.ptr(len, sum);

                for (int d = 0; d <= 9; ++d) {
                    const int nsum = sum + d;
                    if (nsum > kMaxSum) {
                        break;
                    }
                    const int shift = (d * pow10mod_[static_cast<std::size_t>(len)]) % kMod;
                    u64* dst = dp_.ptr(len + 1, nsum);
                    rotate_or(src, dst, shift);
                    dp_.set_nonempty(len + 1, nsum);
                }
            }
        }
    }

    static Node node_from_y(u64 y) {
        Node node;
        node.valid = true;

        node.cnt[9] = y / static_cast<u64>(kMod);
        u64 rem = y % static_cast<u64>(kMod);

        int fac = kMod;
        for (int d = 8; d >= 1; --d) {
            fac /= (d + 1);
            node.cnt[static_cast<std::size_t>(d)] = rem / static_cast<u64>(fac);
            rem %= static_cast<u64>(fac);
        }

        return node;
    }

    DpTable dp_;
    std::array<int, kMaxLen + 1> pow10mod_{};
};

u64 to_small_number(const Node& node) {
    const u64 len = length_of_node(node);
    if (len > 18) {
        return std::numeric_limits<u64>::max();
    }

    u64 value = 0;
    for (int d = 1; d <= 9; ++d) {
        for (u64 c = 0; c < node.cnt[static_cast<std::size_t>(d)]; ++c) {
            value = value * 10ULL + static_cast<u64>(d);
        }
    }
    return value;
}

bool run_checkpoints() {
    Solver solver;

    const Node g5 = solver.best_node_for_sf(5);
    if (to_small_number(g5) != 25ULL) {
        std::cerr << "Checkpoint failed: g(5) should be 25" << '\n';
        return false;
    }

    const Node g20 = solver.best_node_for_sf(20);
    if (to_small_number(g20) != 267ULL) {
        std::cerr << "Checkpoint failed: g(20) should be 267" << '\n';
        return false;
    }

    u64 sum_20 = 0;
    for (int i = 1; i <= 20; ++i) {
        sum_20 += digit_sum_of_node(solver.best_node_for_sf(i));
    }
    if (sum_20 != 156ULL) {
        std::cerr << "Checkpoint failed: sum_{i=1..20} sg(i) should be 156, got " << sum_20 << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    Solver solver;
    std::cout << solver.solve_sum_sg(options.max_i) << '\n';
    return 0;
}
