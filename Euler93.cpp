#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    bool run_checkpoints = true;
};

struct Fraction {
    i64 n = 0;
    i64 d = 1;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

i64 gcd(i64 a, i64 b) {
    while (b != 0) {
        const i64 t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

Fraction reduce(Fraction f) {
    if (f.d < 0) {
        f.n = -f.n;
        f.d = -f.d;
    }
    const i64 g = gcd(f.n, f.d);
    f.n /= g;
    f.d /= g;
    return f;
}

bool apply_op(const Fraction& a, const Fraction& b, const char op, Fraction& out) {
    if (op == '+') {
        out = reduce({a.n * b.d + b.n * a.d, a.d * b.d});
        return true;
    }
    if (op == '-') {
        out = reduce({a.n * b.d - b.n * a.d, a.d * b.d});
        return true;
    }
    if (op == '*') {
        out = reduce({a.n * b.n, a.d * b.d});
        return true;
    }
    if (op == '/') {
        if (b.n == 0) {
            return false;
        }
        out = reduce({a.n * b.d, a.d * b.n});
        return true;
    }
    return false;
}

int consecutive_length_for_digits(const std::array<int, 4>& digits) {
    const std::array<char, 4> ops = {'+', '-', '*', '/'};
    std::set<int> results;

    std::array<int, 4> perm = digits;
    std::sort(perm.begin(), perm.end());

    do {
        const Fraction a{perm[0], 1};
        const Fraction b{perm[1], 1};
        const Fraction c{perm[2], 1};
        const Fraction d{perm[3], 1};

        for (char o1 : ops) {
            for (char o2 : ops) {
                for (char o3 : ops) {
                    Fraction x, y, z;

                    // ((a o1 b) o2 c) o3 d
                    if (apply_op(a, b, o1, x) && apply_op(x, c, o2, y) && apply_op(y, d, o3, z) && z.d == 1 && z.n > 0) {
                        results.insert(static_cast<int>(z.n));
                    }
                    // (a o1 (b o2 c)) o3 d
                    if (apply_op(b, c, o2, x) && apply_op(a, x, o1, y) && apply_op(y, d, o3, z) && z.d == 1 && z.n > 0) {
                        results.insert(static_cast<int>(z.n));
                    }
                    // a o1 ((b o2 c) o3 d)
                    if (apply_op(b, c, o2, x) && apply_op(x, d, o3, y) && apply_op(a, y, o1, z) && z.d == 1 && z.n > 0) {
                        results.insert(static_cast<int>(z.n));
                    }
                    // a o1 (b o2 (c o3 d))
                    if (apply_op(c, d, o3, x) && apply_op(b, x, o2, y) && apply_op(a, y, o1, z) && z.d == 1 && z.n > 0) {
                        results.insert(static_cast<int>(z.n));
                    }
                    // (a o1 b) o2 (c o3 d)
                    if (apply_op(a, b, o1, x) && apply_op(c, d, o3, y) && apply_op(x, y, o2, z) && z.d == 1 && z.n > 0) {
                        results.insert(static_cast<int>(z.n));
                    }
                }
            }
        }
    } while (std::next_permutation(perm.begin(), perm.end()));

    int length = 1;
    while (results.count(length) > 0) {
        ++length;
    }
    return length - 1;
}

std::string solve() {
    int best_len = -1;
    std::array<int, 4> best_digits = {1, 2, 3, 4};

    for (int a = 1; a <= 9; ++a) {
        for (int b = a + 1; b <= 9; ++b) {
            for (int c = b + 1; c <= 9; ++c) {
                for (int d = c + 1; d <= 9; ++d) {
                    const std::array<int, 4> digits = {a, b, c, d};
                    const int len = consecutive_length_for_digits(digits);
                    if (len > best_len) {
                        best_len = len;
                        best_digits = digits;
                    }
                }
            }
        }
    }

    return std::to_string(best_digits[0]) + std::to_string(best_digits[1]) + std::to_string(best_digits[2]) +
           std::to_string(best_digits[3]);
}

bool run_checkpoints() {
    const std::array<int, 4> d = {1, 2, 3, 4};
    if (consecutive_length_for_digits(d) != 28) {
        std::cerr << "Checkpoint failed for digits 1,2,3,4" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
