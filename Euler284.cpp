#include <array>
#include <cstdint>
#include <iostream>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

constexpr int kBase = 14;
constexpr int kDefaultN = 10000;

struct Options {
    int n_max = kDefaultN;
    bool run_checkpoints = true;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n-max=", options.n_max)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n_max >= 1;
}

int mod14(const cpp_int& v) {
    cpp_int r = v % kBase;
    if (r < 0) {
        r += kBase;
    }
    return static_cast<int>(r);
}

std::string to_base14(u64 value) {
    constexpr std::array<char, 14> digits = {
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd'};
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        out.push_back(digits[static_cast<std::size_t>(value % 14)]);
        value /= 14;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

struct SequenceState {
    cpp_int x;
    cpp_int f;
    cpp_int mod;
    int leading_digit;
    u64 digit_sum;

    void step(const std::array<int, kBase>& inv_mod14) {
        const cpp_int x_old = x;
        const cpp_int mod_old = mod;

        const int a = mod14(2 * x_old - 1);
        const int inv = inv_mod14[static_cast<std::size_t>(a)];
        const int f_mod = mod14(f);
        const int t = (kBase - f_mod) % kBase * inv % kBase;

        x = x_old + static_cast<cpp_int>(t) * mod_old;
        f = (f + (2 * x_old - 1) * t + static_cast<cpp_int>(t) * t * mod_old) / kBase;
        mod = mod_old * kBase;

        leading_digit = t;
        digit_sum += static_cast<u64>(t);
    }
};

SequenceState make_initial_state(int root) {
    SequenceState s;
    s.x = root;
    s.mod = kBase;
    s.f = (s.x * s.x - s.x) / kBase;
    s.leading_digit = root;
    s.digit_sum = static_cast<u64>(root);
    return s;
}

u64 solve_sum_of_digit_sums(const int n_max) {
    std::array<int, kBase> inv_mod14{};
    inv_mod14.fill(-1);
    for (int a = 1; a < kBase; ++a) {
        for (int b = 1; b < kBase; ++b) {
            if ((a * b) % kBase == 1) {
                inv_mod14[static_cast<std::size_t>(a)] = b;
                break;
            }
        }
    }

    SequenceState seq7 = make_initial_state(7);
    SequenceState seq8 = make_initial_state(8);

    u64 total = 1;  // steady square "1" for n = 1

    for (int n = 1; n <= n_max; ++n) {
        if (seq7.leading_digit != 0) {
            total += seq7.digit_sum;
        }
        if (seq8.leading_digit != 0) {
            total += seq8.digit_sum;
        }

        if (n == n_max) {
            break;
        }

        seq7.step(inv_mod14);
        seq8.step(inv_mod14);
    }

    return total;
}

bool run_checkpoints() {
    if (solve_sum_of_digit_sums(9) != 582ULL) {
        std::cerr << "Checkpoint failed for n<=9" << '\n';
        return false;
    }

    // Verify steady-square condition during the first few lifts.
    SequenceState seq7 = make_initial_state(7);
    SequenceState seq8 = make_initial_state(8);
    std::array<int, kBase> inv_mod14{};
    inv_mod14.fill(-1);
    for (int a = 1; a < kBase; ++a) {
        for (int b = 1; b < kBase; ++b) {
            if ((a * b) % kBase == 1) {
                inv_mod14[static_cast<std::size_t>(a)] = b;
                break;
            }
        }
    }

    for (int n = 1; n <= 20; ++n) {
        if ((seq7.x * seq7.x - seq7.x) % seq7.mod != 0) {
            std::cerr << "Idempotence checkpoint failed for root 7 at n=" << n << '\n';
            return false;
        }
        if ((seq8.x * seq8.x - seq8.x) % seq8.mod != 0) {
            std::cerr << "Idempotence checkpoint failed for root 8 at n=" << n << '\n';
            return false;
        }
        seq7.step(inv_mod14);
        seq8.step(inv_mod14);
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

    const u64 answer = solve_sum_of_digit_sums(options.n_max);
    std::cout << to_base14(answer) << '\n';
    return 0;
}
