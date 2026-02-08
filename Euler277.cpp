#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = long long;
using i128 = __int128_t;

struct Options {
    std::string sequence = "UDDDUdddDDUDDddDdDddDDUDDdUUDd";
    i64 lower_bound = 1000000000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_string_after_prefix(const std::string& arg, const std::string& prefix, std::string& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    value = arg.substr(prefix.size());
    return !value.empty();
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_i64_after_prefix(arg, "--lower-bound=", options.lower_bound) ||
            parse_string_after_prefix(arg, "--sequence=", options.sequence)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.lower_bound < 0) {
        return false;
    }
    for (char c : options.sequence) {
        if (c != 'U' && c != 'D' && c != 'd') {
            return false;
        }
    }
    return true;
}

i64 gcd_i64(i64 a, i64 b) {
    while (b != 0) {
        const i64 t = a % b;
        a = b;
        b = t;
    }
    return a >= 0 ? a : -a;
}

i64 mod_inverse(i64 a, i64 mod) {
    i64 t = 0;
    i64 new_t = 1;
    i64 r = mod;
    i64 new_r = a;

    while (new_r != 0) {
        const i64 q = r / new_r;
        const i64 next_t = t - q * new_t;
        t = new_t;
        new_t = next_t;

        const i64 next_r = r - q * new_r;
        r = new_r;
        new_r = next_r;
    }

    if (r != 1) {
        return -1;
    }

    t %= mod;
    if (t < 0) {
        t += mod;
    }
    return t;
}

i128 solve(const std::string& sequence, const i128 lower_bound) {
    i128 A = 1;
    i128 B = 0;
    i128 C = 1;

    for (int i = static_cast<int>(sequence.size()) - 1; i >= 0; --i) {
        const char step = sequence[static_cast<std::size_t>(i)];
        if (step == 'D') {
            A = 3 * A;
            B = 3 * B;
        } else if (step == 'U') {
            A = 3 * A;
            B = 3 * B - 2 * C;
            C = 4 * C;
        } else {
            A = 3 * A;
            B = 3 * B + C;
            C = 2 * C;
        }
    }

    const i64 a64 = static_cast<i64>(A);
    const i64 c64 = static_cast<i64>(C);
    const i64 g = gcd_i64(a64, c64);

    const i64 mod = c64 / g;
    const i64 a_reduced = a64 / g;

    const i64 a_mod = (a_reduced % mod + mod) % mod;
    const i64 inv = mod_inverse(a_mod, mod);
    if (inv < 0) {
        return -1;
    }

    const i128 B_reduced = B / g;
    const i64 rhs = static_cast<i64>(((-B_reduced) % mod + mod) % mod);
    const i64 t0 = static_cast<i64>((static_cast<i128>(rhs) * inv) % mod);

    const i128 n0 = (A * static_cast<i128>(t0) + B) / C;
    const i128 n_step = A / g;
    const i128 t_step = mod;

    i128 k = 0;
    if (n0 <= lower_bound) {
        k = (lower_bound - n0) / n_step + 1;
    }

    const i128 t_min = static_cast<i128>(t0) + k * t_step;
    if (t_min <= 0) {
        const i128 add = (-t_min) / t_step + 1;
        k += add;
    }

    const i128 t = static_cast<i128>(t0) + k * t_step;
    return (A * t + B) / C;
}

bool follows_prefix(i128 value, const std::string& sequence) {
    for (char c : sequence) {
        const int mod3 = static_cast<int>(value % 3);
        if (c == 'D') {
            if (mod3 != 0) {
                return false;
            }
            value /= 3;
        } else if (c == 'U') {
            if (mod3 != 1) {
                return false;
            }
            value = (4 * value + 2) / 3;
        } else {
            if (mod3 != 2) {
                return false;
            }
            value = (2 * value - 1) / 3;
        }
    }
    return true;
}

bool run_checkpoints() {
    const std::string sample_seq = "DdDddUUdDD";
    const i128 sample = solve(sample_seq, 1000000);
    if (sample != 1004064) {
        std::cerr << "Checkpoint failed for sample sequence" << '\n';
        return false;
    }
    if (!follows_prefix(sample, sample_seq)) {
        std::cerr << "Checkpoint failed for sample prefix validation" << '\n';
        return false;
    }
    return true;
}

std::string i128_to_string(i128 value) {
    if (value == 0) {
        return "0";
    }
    bool neg = value < 0;
    if (neg) {
        value = -value;
    }

    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (neg) {
        out.push_back('-');
    }
    std::reverse(out.begin(), out.end());
    return out;
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

    const i128 answer = solve(options.sequence, options.lower_bound);
    std::cout << i128_to_string(answer) << '\n';
    return 0;
}
