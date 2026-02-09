#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int k_max = 13;
    bool run_checkpoints = true;
    u64 single_n = 0;
};

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

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<u64>(c - '0');
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
        if (parse_int_after_prefix(arg, "--k-max=", options.k_max)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--single=", options.single_n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.k_max >= 1;
}

struct Op {
    std::vector<std::uint8_t> end_state;
    std::vector<u64> add_count;

    bool operator==(const Op& other) const {
        return end_state == other.end_state && add_count == other.add_count;
    }
};

struct OpHash {
    std::size_t operator()(const Op& op) const {
        std::size_t h = 1469598103934665603ULL;
        for (std::uint8_t v : op.end_state) {
            h ^= static_cast<std::size_t>(v) + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
        }
        for (u64 v : op.add_count) {
            h ^= static_cast<std::size_t>(v ^ (v >> 32U)) + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
        }
        return h;
    }
};

class PatternEngine {
public:
    explicit PatternEngine(const std::string& pattern)
        : pattern_(pattern), L_(static_cast<int>(pattern.size())) {
        build_automaton();
        build_digit_operators();
        build_pow10(20);
    }

    u64 nth_occurrence_position(const u64 target_occurrence) {
        u64 hi = 1;
        while (occ_and_state(hi).first < target_occurrence) {
            if (hi > std::numeric_limits<u64>::max() / 2ULL) {
                hi = std::numeric_limits<u64>::max();
                break;
            }
            hi *= 2ULL;
        }

        u64 lo = 1;
        while (lo < hi) {
            const u64 mid = lo + (hi - lo) / 2ULL;
            if (occ_and_state(mid).first >= target_occurrence) {
                hi = mid;
            } else {
                lo = mid + 1ULL;
            }
        }

        const u64 number = lo;
        const auto [before_count, before_state] = occ_and_state(number - 1ULL);
        const u64 need_inside = target_occurrence - before_count;

        u64 prefix_digits = digits_before(number - 1ULL);
        u64 found_inside = 0;
        int state = before_state;

        const std::string num = std::to_string(number);
        for (std::size_t i = 0; i < num.size(); ++i) {
            const int d = num[i] - '0';
            if (out_[state][d] != 0) {
                ++found_inside;
                if (found_inside == need_inside) {
                    const u64 end_pos = prefix_digits + static_cast<u64>(i + 1);
                    return end_pos - static_cast<u64>(L_) + 1ULL;
                }
            }
            state = next_[state][d];
        }

        return 0;
    }

private:
    std::string pattern_;
    int L_ = 0;

    std::vector<std::vector<int>> next_;
    std::vector<std::vector<std::uint8_t>> out_;

    std::vector<Op> ops_;
    std::unordered_map<Op, int, OpHash> op_to_id_;

    int id_identity_ = -1;
    std::array<int, 10> digit_op_id_{};

    std::unordered_map<u64, int> compose_cache_;
    std::unordered_map<u64, int> append_cache_;
    std::unordered_map<u64, int> full_block_cache_;

    std::vector<u64> pow10_;

    static u64 make_key(const int a, const int b) {
        return (static_cast<u64>(a) << 32U) | static_cast<std::uint32_t>(b);
    }

    int intern_op(const Op& op) {
        const auto it = op_to_id_.find(op);
        if (it != op_to_id_.end()) {
            return it->second;
        }
        const int id = static_cast<int>(ops_.size());
        ops_.push_back(op);
        op_to_id_.emplace(ops_[static_cast<std::size_t>(id)], id);
        return id;
    }

    void build_automaton() {
        std::vector<int> pi(static_cast<std::size_t>(L_), 0);
        for (int i = 1; i < L_; ++i) {
            int j = pi[static_cast<std::size_t>(i - 1)];
            while (j > 0 && pattern_[static_cast<std::size_t>(i)] != pattern_[static_cast<std::size_t>(j)]) {
                j = pi[static_cast<std::size_t>(j - 1)];
            }
            if (pattern_[static_cast<std::size_t>(i)] == pattern_[static_cast<std::size_t>(j)]) {
                ++j;
            }
            pi[static_cast<std::size_t>(i)] = j;
        }

        next_.assign(static_cast<std::size_t>(L_ + 1), std::vector<int>(10, 0));
        out_.assign(static_cast<std::size_t>(L_ + 1), std::vector<std::uint8_t>(10, 0));

        for (int st = 0; st <= L_; ++st) {
            for (int d = 0; d <= 9; ++d) {
                const char c = static_cast<char>('0' + d);
                int j = st;
                while (j > 0 && (j == L_ || pattern_[static_cast<std::size_t>(j)] != c)) {
                    j = pi[static_cast<std::size_t>(j - 1)];
                }
                if (j < L_ && pattern_[static_cast<std::size_t>(j)] == c) {
                    ++j;
                }

                if (j == L_) {
                    out_[static_cast<std::size_t>(st)][static_cast<std::size_t>(d)] = 1;
                    j = pi[static_cast<std::size_t>(L_ - 1)];
                }
                next_[static_cast<std::size_t>(st)][static_cast<std::size_t>(d)] = j;
            }
        }
    }

    void build_digit_operators() {
        Op identity;
        identity.end_state.resize(static_cast<std::size_t>(L_ + 1));
        identity.add_count.assign(static_cast<std::size_t>(L_ + 1), 0);
        for (int s = 0; s <= L_; ++s) {
            identity.end_state[static_cast<std::size_t>(s)] = static_cast<std::uint8_t>(s);
        }
        id_identity_ = intern_op(identity);

        for (int d = 0; d <= 9; ++d) {
            Op op;
            op.end_state.resize(static_cast<std::size_t>(L_ + 1));
            op.add_count.resize(static_cast<std::size_t>(L_ + 1));
            for (int s = 0; s <= L_; ++s) {
                op.end_state[static_cast<std::size_t>(s)] =
                    static_cast<std::uint8_t>(next_[static_cast<std::size_t>(s)][static_cast<std::size_t>(d)]);
                op.add_count[static_cast<std::size_t>(s)] =
                    static_cast<u64>(out_[static_cast<std::size_t>(s)][static_cast<std::size_t>(d)]);
            }
            digit_op_id_[static_cast<std::size_t>(d)] = intern_op(op);
        }
    }

    void build_pow10(const int max_power) {
        pow10_.assign(static_cast<std::size_t>(max_power + 1), 1ULL);
        for (int i = 1; i <= max_power; ++i) {
            pow10_[static_cast<std::size_t>(i)] = pow10_[static_cast<std::size_t>(i - 1)] * 10ULL;
        }
    }

    int compose_ids(const int id_a, const int id_b) {
        const u64 key = make_key(id_a, id_b);
        const auto it = compose_cache_.find(key);
        if (it != compose_cache_.end()) {
            return it->second;
        }

        const Op& A = ops_[static_cast<std::size_t>(id_a)];
        const Op& B = ops_[static_cast<std::size_t>(id_b)];

        Op C;
        C.end_state.resize(static_cast<std::size_t>(L_ + 1));
        C.add_count.resize(static_cast<std::size_t>(L_ + 1));

        for (int s = 0; s <= L_; ++s) {
            const int mid = A.end_state[static_cast<std::size_t>(s)];
            C.end_state[static_cast<std::size_t>(s)] = B.end_state[static_cast<std::size_t>(mid)];
            C.add_count[static_cast<std::size_t>(s)] =
                A.add_count[static_cast<std::size_t>(s)] + B.add_count[static_cast<std::size_t>(mid)];
        }

        const int id = intern_op(C);
        compose_cache_.emplace(key, id);
        return id;
    }

    int append_digit(const int op_id, const int d) {
        const u64 key = make_key(op_id, d);
        const auto it = append_cache_.find(key);
        if (it != append_cache_.end()) {
            return it->second;
        }
        const int id = compose_ids(op_id, digit_op_id_[static_cast<std::size_t>(d)]);
        append_cache_.emplace(key, id);
        return id;
    }

    int full_block(const int prefix_op, const int rem_digits) {
        const u64 key = make_key(prefix_op, rem_digits);
        const auto it = full_block_cache_.find(key);
        if (it != full_block_cache_.end()) {
            return it->second;
        }

        int id;
        if (rem_digits == 0) {
            id = prefix_op;
        } else {
            int result = id_identity_;
            for (int d = 0; d <= 9; ++d) {
                const int child = full_block(append_digit(prefix_op, d), rem_digits - 1);
                result = compose_ids(result, child);
            }
            id = result;
        }

        full_block_cache_.emplace(key, id);
        return id;
    }

    int partial_block(const int prefix_op, const int rem_digits, const u64 upper_suffix) {
        if (rem_digits == 0) {
            return prefix_op;
        }

        const u64 base = pow10_[static_cast<std::size_t>(rem_digits - 1)];
        const int first = static_cast<int>(upper_suffix / base);
        const u64 rest = upper_suffix % base;

        int result = id_identity_;
        for (int d = 0; d < first; ++d) {
            const int child = full_block(append_digit(prefix_op, d), rem_digits - 1);
            result = compose_ids(result, child);
        }

        const int last = partial_block(append_digit(prefix_op, first), rem_digits - 1, rest);
        result = compose_ids(result, last);
        return result;
    }

    int range_operator(const u64 N) {
        if (N == 0) {
            return id_identity_;
        }

        const std::string s = std::to_string(N);
        const int D = static_cast<int>(s.size());

        int result = id_identity_;

        for (int d = 1; d < D; ++d) {
            const int rem = d - 1;
            for (int a = 1; a <= 9; ++a) {
                const int block = full_block(digit_op_id_[static_cast<std::size_t>(a)], rem);
                result = compose_ids(result, block);
            }
        }

        const int first = s[0] - '0';
        const int rem = D - 1;

        for (int a = 1; a < first; ++a) {
            const int block = full_block(digit_op_id_[static_cast<std::size_t>(a)], rem);
            result = compose_ids(result, block);
        }

        if (rem == 0) {
            result = compose_ids(result, digit_op_id_[static_cast<std::size_t>(first)]);
        } else {
            u64 suffix = 0;
            for (int i = 1; i < D; ++i) {
                suffix = suffix * 10ULL + static_cast<u64>(s[static_cast<std::size_t>(i)] - '0');
            }
            const int part = partial_block(digit_op_id_[static_cast<std::size_t>(first)], rem, suffix);
            result = compose_ids(result, part);
        }

        return result;
    }

    std::pair<u64, int> occ_and_state(const u64 N) {
        const int op_id = range_operator(N);
        const Op& op = ops_[static_cast<std::size_t>(op_id)];
        return {op.add_count[0], static_cast<int>(op.end_state[0])};
    }

    u64 digits_before(const u64 N) const {
        if (N == 0) {
            return 0;
        }

        const std::string s = std::to_string(N);
        const int D = static_cast<int>(s.size());

        u64 total = 0;
        for (int d = 1; d < D; ++d) {
            total += 9ULL * pow10_[static_cast<std::size_t>(d - 1)] * static_cast<u64>(d);
        }

        total += (N - pow10_[static_cast<std::size_t>(D - 1)] + 1ULL) * static_cast<u64>(D);
        return total;
    }
};

u64 brute_nth_position_small(const int n) {
    const std::string token = std::to_string(n);
    const int L = static_cast<int>(token.size());

    std::string stream;
    stream.reserve(10000);
    std::vector<u64> starts;
    starts.reserve(static_cast<std::size_t>(n));

    u64 checked_until = 0;

    for (int x = 1; static_cast<int>(starts.size()) < n; ++x) {
        const std::string s = std::to_string(x);
        const u64 old_size = static_cast<u64>(stream.size());
        stream += s;

        const u64 begin = (old_size >= static_cast<u64>(L - 1)) ? (old_size - static_cast<u64>(L - 1)) : 0ULL;
        for (u64 i = begin; i + static_cast<u64>(L) <= static_cast<u64>(stream.size()); ++i) {
            bool match = true;
            for (int j = 0; j < L; ++j) {
                if (stream[static_cast<std::size_t>(i + static_cast<u64>(j))] != token[static_cast<std::size_t>(j)]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                if (starts.empty() || starts.back() != i + 1ULL) {
                    starts.push_back(i + 1ULL);
                }
            }
        }
        checked_until = old_size;
        (void)checked_until;
    }

    return starts[static_cast<std::size_t>(n - 1)];
}

bool run_checkpoints() {
    struct Sample {
        u64 n;
        u64 expected;
    };

    const std::vector<Sample> samples = {
        {1ULL, 1ULL},
        {5ULL, 81ULL},
        {12ULL, 271ULL},
        {7780ULL, 111111365ULL},
    };

    for (const auto& s : samples) {
        PatternEngine engine(std::to_string(s.n));
        const u64 got = engine.nth_occurrence_position(s.n);
        if (got != s.expected) {
            std::cerr << "Sample checkpoint failed for n=" << s.n << ": got " << got
                      << ", expected " << s.expected << '\n';
            return false;
        }
    }

    for (int n = 1; n <= 30; ++n) {
        PatternEngine engine(std::to_string(n));
        const u64 fast = engine.nth_occurrence_position(static_cast<u64>(n));
        const u64 brute = brute_nth_position_small(n);
        if (fast != brute) {
            std::cerr << "Brute checkpoint failed for n=" << n << ": fast=" << fast
                      << ", brute=" << brute << '\n';
            return false;
        }
    }

    return true;
}

u64 solve_sum_k(const int k_max) {
    u64 sum = 0;
    u64 n = 1;
    for (int k = 1; k <= k_max; ++k) {
        n *= 3ULL;
        PatternEngine engine(std::to_string(n));
        sum += engine.nth_occurrence_position(n);
    }
    return sum;
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

    if (options.single_n > 0) {
        PatternEngine engine(std::to_string(options.single_n));
        std::cout << engine.nth_occurrence_position(options.single_n) << '\n';
        return 0;
    }

    std::cout << solve_sum_k(options.k_max) << '\n';
    return 0;
}
