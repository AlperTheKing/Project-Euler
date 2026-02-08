#include <algorithm>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr int kDefaultM = 5000000;
constexpr int kExponentLimit = 12000;

struct Options {
    int m = kDefaultM;
    bool run_checkpoints = true;
};

bool parse_nonnegative_int(const std::string& text, int& out) {
    if (text.empty()) {
        return false;
    }

    std::uint64_t value = 0;
    for (const char ch : text) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        value = value * 10ULL + static_cast<std::uint64_t>(ch - '0');
        if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    out = static_cast<int>(value);
    return true;
}

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    return parse_nonnegative_int(arg.substr(prefix.size()), value);
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    bool seen_positional_m = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        int parsed = 0;
        if (parse_int_after_prefix(arg, "--m=", parsed)) {
            options.m = parsed;
            continue;
        }

        if (!seen_positional_m && parse_nonnegative_int(arg, parsed)) {
            options.m = parsed;
            seen_positional_m = true;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.m <= 0) {
        std::cerr << "m must be positive.\n";
        return false;
    }

    return true;
}

class OddRepresentation {
public:
    explicit OddRepresentation(const int exponent_limit)
        : offset_(exponent_limit),
          bit_(2 * exponent_limit + 1, 0),
          position_(2 * exponent_limit + 1, -1) {}

    int offset() const { return offset_; }
    int size() const { return static_cast<int>(bit_.size()); }
    int weight() const { return static_cast<int>(active_indices_.size()); }

    bool has(const int index) const {
        if (index < 0 || index >= static_cast<int>(bit_.size())) {
            return false;
        }
        return bit_[index] != 0;
    }

    void set(const int index) {
        if (bit_[index] != 0) {
            return;
        }
        bit_[index] = 1;
        position_[index] = static_cast<int>(active_indices_.size());
        active_indices_.push_back(index);
    }

    void clear(const int index) {
        if (index < 0 || index >= static_cast<int>(bit_.size()) ||
            bit_[index] == 0) {
            return;
        }

        bit_[index] = 0;
        const int pos = position_[index];
        const int last_index = active_indices_.back();
        active_indices_[pos] = last_index;
        position_[last_index] = pos;
        active_indices_.pop_back();
        position_[index] = -1;
    }

    const std::vector<int>& active_indices() const { return active_indices_; }

private:
    int offset_;
    std::vector<std::uint8_t> bit_;
    std::vector<int> position_;
    std::vector<int> active_indices_;
};

class SquareRepresentation {
public:
    explicit SquareRepresentation(const int exponent_limit)
        : offset_(exponent_limit), bit_(2 * exponent_limit + 1, 0) {}

    int offset() const { return offset_; }
    int size() const { return static_cast<int>(bit_.size()); }
    int weight() const { return weight_; }

    bool has(const int index) const {
        if (index < 0 || index >= static_cast<int>(bit_.size())) {
            return false;
        }
        return bit_[index] != 0;
    }

    void set(const int index) {
        if (bit_[index] != 0) {
            return;
        }
        bit_[index] = 1;
        ++weight_;
    }

    void clear(const int index) {
        if (index < 0 || index >= static_cast<int>(bit_.size()) ||
            bit_[index] == 0) {
            return;
        }
        bit_[index] = 0;
        --weight_;
    }

private:
    int offset_;
    std::vector<std::uint8_t> bit_;
    int weight_ = 0;
};

template <typename Representation>
void add_single_power(Representation& rep, const int index) {
    static thread_local std::vector<int> stack;
    stack.clear();
    stack.push_back(index);

    const int n = rep.size();
    while (!stack.empty()) {
        const int x = stack.back();
        stack.pop_back();

        // Rules below may access up to x+3 and x-7.
        if (x < 8 || x + 8 >= n) {
            throw std::runtime_error(
                "Exponent window exceeded; increase kExponentLimit.");
        }

        if (rep.has(x + 2)) {
            rep.clear(x + 2);
            stack.push_back(x + 3);
            continue;
        }
        if (rep.has(x - 2)) {
            rep.clear(x - 2);
            stack.push_back(x + 1);
            continue;
        }
        if (rep.has(x - 1)) {
            rep.clear(x - 1);
            stack.push_back(x - 4);
            stack.push_back(x + 1);
            continue;
        }
        if (rep.has(x + 1)) {
            rep.clear(x + 1);
            stack.push_back(x - 3);
            stack.push_back(x + 2);
            continue;
        }
        if (rep.has(x)) {
            rep.clear(x);
            stack.push_back(x - 7);
            stack.push_back(x - 2);
            stack.push_back(x + 1);
            continue;
        }

        rep.set(x);
    }
}

std::vector<int> sorted_exponents(const OddRepresentation& rep) {
    std::vector<int> result;
    result.reserve(rep.active_indices().size());
    for (const int idx : rep.active_indices()) {
        result.push_back(idx - rep.offset());
    }
    std::sort(result.begin(), result.end());
    return result;
}

bool validate_integer_examples() {
    OddRepresentation value(kExponentLimit);

    std::vector<int> exponents_for_3;
    std::vector<int> exponents_for_10;

    for (int n = 1; n <= 10; ++n) {
        add_single_power(value, value.offset());

        if (n == 3) {
            if (value.weight() != 4) {
                std::cerr << "Validation failed: w(3) should be 4, got "
                          << value.weight() << "\n";
                return false;
            }
            exponents_for_3 = sorted_exponents(value);
        }

        if (n == 10) {
            if (value.weight() != 3) {
                std::cerr << "Validation failed: w(10) should be 3, got "
                          << value.weight() << "\n";
                return false;
            }
            exponents_for_10 = sorted_exponents(value);
        }
    }

    const std::vector<int> expected_3 = {-10, -5, -1, 2};
    const std::vector<int> expected_10 = {-10, -7, 6};
    if (exponents_for_3 != expected_3) {
        std::cerr << "Validation failed: exponent set for 3 mismatch.\n";
        return false;
    }
    if (exponents_for_10 != expected_10) {
        std::cerr << "Validation failed: exponent set for 10 mismatch.\n";
        return false;
    }

    return true;
}

struct SolveResult {
    u64 s_m = 0;
    u64 s_10 = 0;
    u64 s_1000 = 0;
};

SolveResult solve_squares(const int m) {
    OddRepresentation odd(kExponentLimit);
    SquareRepresentation square(kExponentLimit);

    SolveResult result;

    for (int i = 1; i <= m; ++i) {
        // odd_i = 2*i - 1. This performs odd += 2 with the first step as +1.
        add_single_power(odd, odd.offset());
        if (i != 1) {
            add_single_power(odd, odd.offset());
        }

        // square_i = square_{i-1} + odd_i.
        for (const int idx : odd.active_indices()) {
            add_single_power(square, idx);
        }

        result.s_m += static_cast<u64>(square.weight());
        if (i == 10) {
            result.s_10 = result.s_m;
        }
        if (i == 1000) {
            result.s_1000 = result.s_m;
        }
    }

    return result;
}

bool run_validations() {
    if (!validate_integer_examples()) {
        return false;
    }

    const SolveResult sample = solve_squares(1000);
    if (sample.s_10 != 61ULL) {
        std::cerr << "Validation failed: S(10) should be 61, got "
                  << sample.s_10 << "\n";
        return false;
    }
    if (sample.s_1000 != 19403ULL) {
        std::cerr << "Validation failed: S(1000) should be 19403, got "
                  << sample.s_1000 << "\n";
        return false;
    }

    return true;
}

}  // namespace

int main(const int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    try {
        if (options.run_checkpoints && !run_validations()) {
            return 1;
        }

        const SolveResult result = solve_squares(options.m);
        std::cout << result.s_m << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
